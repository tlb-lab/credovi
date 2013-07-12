import os
from itertools import groupby
from operator import itemgetter
from collections import OrderedDict

from sqlalchemy import text
from openeye.oechem import oemolistream, OECount, OEFormat_OEB, OEIsHeavy
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

from credovi import app
from credovi.schema import engine
from credovi.util.timer import Timer
from credovi.structbio import structure as struct
from eyesopen.oechem import mol_to_smiles, mol_to_oeb, mol_to_pdb, mol_to_sdf
from credovi.schema.tables.ligands import ligand_molstrings, ligand_usr

def get_ligand_ids():
    """
    Returns a mapping between new ligand identifiers and their PDB identifiers.
    This mapping is necessary to link ligand identifiers to the ligands that are
    attached to structures.
    """
    statement = text("""
                   SELECT s.pdb, b.assembly_serial, l.entity_serial, l.ligand_id
                     FROM {schema}.ligands l
                     JOIN {schema}.biomolecules b USING(biomolecule_id)
                     JOIN {schema}.structures s USING(structure_id)
                LEFT JOIN {schema}.ligand_molstrings lm USING(ligand_id)
                    WHERE lm.ligand_id IS NULL
                          AND l.num_hvy_atoms >= :min_hvy_atoms
                 ORDER BY 1,2,3,4
                """.format(schema=app.config.get('database','schema')))

    result = engine.execute(statement, min_hvy_atoms=5).fetchall()

    mapping = OrderedDict()

    for pdb, biomoliter in groupby(result, key=itemgetter(0)):
        mapping[pdb] = {}

        for biomolecule, groupiter in groupby(biomoliter, key=itemgetter(1)):
            mapping[pdb].update({biomolecule:{}})

            for pdb, biomolecule, entity_id, ligand_id in groupiter:
                mapping[pdb][biomolecule].update({entity_id:ligand_id})

    return mapping

def do(controller):
    """
    """
    # timer to clock functions and parts of the program
    timer = Timer()
    timer.start("app")

    # get the controller command
    cmd = controller.command

    # get the command line arguments and options
    args = controller.pargs

    insert = ligand_molstrings.insert()

    ifs = oemolistream()
    ifs.SetFormat(OEFormat_OEB)

    # directory containing all the biological assemblies in OEB format
    QUAT_OEB_DIR = app.config.get('directories','quat_oeb')

    data = get_ligand_ids()

    app.log.debug("retrieved {} PDB entries from CREDO without ligand molstrings."
                  .format(len(data)))

    # initialize progressbar
    if args.progressbar:
        bar = ProgressBar(widgets=['PDB entries: ', SimpleProgress(), ' ',
                                   Percentage(), Bar()],
                          maxval=len(data.keys())).start()

    for counter, pdb in enumerate(data,1):
        if args.progressbar: bar.update(counter)

        for assembly_serial in data[pdb]:

            # get the path to the biological assembly containing the ligands
            filename = "{}-{}.oeb".format(pdb.lower(), assembly_serial)
            path = os.path.join(QUAT_OEB_DIR, pdb[1:3].lower(), pdb.lower(),
                                filename)

            if not os.path.isfile: app.log.warn("file {path} does not exist!"
                                                .format(path=path))

            # load the structure as OEGraphMol
            # returns None if structure could not be loaded
            biomolecule = struct.get_structure(str(path), str(pdb))

            if not biomolecule:
                app.log.warn("cannot read biomolecule {}!".format(path))
                continue

            # biomolecule does not have any attached ligands, this should not happen
            if not biomolecule.HasData('ligands'):
                app.log.warn("cannot insert molstrings for {}: assembly does not "
                             "have any attached ligands!".format(path))
                continue

            # method can return either scalar or list
            ligands = biomolecule.GetListData('ligands')

            for ligand in ligands:
                entity_serial = ligand.GetIntData('entity_serial')

                # ignore ligands that are too small
                if OECount(ligand, OEIsHeavy()) < 5: continue

                # get the corresponding credo ligand identifier for this ligand
                try: ligand_id = data[pdb][assembly_serial][entity_serial]

                # this ligand is not of our interest; either too small or already
                # in the molstrings table
                except KeyError:
                    app.log.warn("cannot insert ligand molstrings for ligand: "
                                 "ligand_id cannot be obtained for {}-{}-{}"
                                 .format(pdb, assembly_serial, entity_serial))
                    continue

                # convert ligand molecule to desired file formats
                ism = mol_to_smiles(ligand, from3d=True, isomeric=True, reset_charges=True)
                pdbformat = mol_to_pdb(ligand)
                oeb = mol_to_oeb(ligand)
                sdf = mol_to_sdf(ligand)

                if not args.dry_run:

                    # insert data into database
                    engine.execute(insert,
                                   ligand_id=ligand_id,
                                   ism=ism, pdb=pdbformat, oeb=oeb, sdf=sdf)

    if args.progressbar: bar.finish()

    # insert the USRCAT moments for the newly inserted ligand structures
    if not args.dry_run and args.usr:

        # insert usr moments for all ligands
        engine.execute("""
                          INSERT INTO {schema}.ligand_usr(ligand_id, usr_space, usr_moments)
                            WITH moments AS
                                 (
                                     SELECT lm.ligand_id, openeye.usrcat(oeb) as moments
                                       FROM {schema}.ligand_molstrings lm
                                       JOIN {schema}.ligands l on l.ligand_id = lm.ligand_id
                                  LEFT JOIN {schema}.ligand_usr lu ON lu.ligand_id = lm.ligand_id
                                      WHERE l.num_hvy_atoms >= 7 AND lu.ligand_id IS NULL
                                 )
                          SELECT ligand_id, cube(moments[1:12]), moments
                            FROM moments
                        ORDER BY 1
                       """.format(schema=app.config.get('database','schema')))

        app.log.debug("inserted USRCAT moments.")

         # insert usr moments for all ligands
        engine.execute("""
                         UPDATE {schema}.ligand_usr lu
                            SET npr1 = (sq.nprs).npr1, npr2 = (sq.nprs).npr2
                           FROM (
                                 SELECT ligand_id, openeye.nprs(oeb) as nprs
                                   FROM {schema}.ligand_molstrings lm
                                ) sq
                          WHERE sq.ligand_id = lu.ligand_id;
                       """.format(schema=app.config.get('database','schema')))

        app.log.debug("Updated NPRs.")

        # cluster table using the gist index
        engine.execute("CLUSTER {schema}.ligand_usr"
                       .format(schema=app.config.get('database','schema')))
