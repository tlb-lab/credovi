"""
This module is used to prepare structures in PDB format for usage with the CREDO
contact generation workflow. Preparation includes the identification of molecular
entities inside the structure, the assignment of atom types flags and the generation
of quaternary assemblies.

python credovi.py mmcif currentpdbs -O 2VQE -L 10000 | parallel --eta --halt 2 -n 4 python credovi.py credo preparepdb -Q --oeb --pdb -S {1},{2},{3},{4}
"""

import os
import sys
import time
import glob
import string

from progressbar import ProgressBar, Percentage, Bar, SimpleProgress


from credovi import app # to get the configuration
from credovi.structbio import db, biomol
from credovi.structbio.structure import (assign_entity_serials, calc_calpha_ratio,
                                          get_ligands, identify_surface_atoms,
                                          parse_header, set_atom_type_flags)
from credovi.lib.openeye import *

# PDB mirror directory
PDB_MIRROR_DIR      = app.config.get('directories','pdb')

# directories where prepared structured will be saved
QUAT_OEB_DIR        = app.config.get('directories','quat_oeb')
QUAT_PDB_DIR        = app.config.get('directories','quat_pdb')
QUAT_OEB_CHAIN_DIR  = app.config.get('directories','quat_chain_oeb')

def get_pdbs_to_process(args):
    """
    Returns a list containing the PDB codes of PDB entries that should be processed.
    """
    pdbs = []
    
    # process all PDB files on mirror
    if args.all:
        pass

    # process only specific pdb entries
    elif args.structures:
        pdbs = [pdb.lower() for pdb in args.structures.strip().split(',') if pdb]

    # use only a limited subset of the pdb
    elif args.testset:
        pdbs = map(string.lower, app.config['test sets'][args.testset])
    
    # apply an offset to the list of PDBs by using the given offset PDB as index
    if pdbs and args.offset:
        try:
            index = pdbs.index(args.offset.lower())
            pdbs = pdbs[index:]
        except ValueError:
            app.log.fatal("cannot apply offset: {0} is not in list of PDBs to "
                          "process.".format(args.offset))
            app.close()
            sys.exit(1)
    
    return sorted(pdbs)

def do(controller):
    """
    """
    # get the command line arguments and options
    args = controller.pargs    
    
    # get the PDB codes of PDB entries that should be processed
    pdbs = get_pdbs_to_process(args)
    
    # OEChem input flavour to use for reading PDB entries
    # this flavour will keep alternate location groups
    OE_PDB_INPUT_FLAVOR = OEIFlavor_Generic_Default | OEIFlavor_PDB_ALTLOC | OEIFlavor_PDB_Default | OEIFlavor_PDB_DATA    
    
    # PDB structure reader for gzipped structures
    ifs = oemolistream()
    ifs.SetFlavor(OEFormat_PDB, OE_PDB_INPUT_FLAVOR)
    ifs.Setgz(True)

    # PDB structure writer
    oebofs = oemolostream()
    oebofs.SetFormat(OEFormat_OEB)

    # ligand structure writer
    pdbofs = oemolostream()
    pdbofs.SetFormat(OEFormat_PDB)

    total, counter = len(pdbs), 0.0

    # initialize progress bar
    if args.progressbar:
        
        # INFO message would disrupt the progress bar
        app.log.set_level('WARN')
        
        bar = ProgressBar(widgets=['PDB entries: ', SimpleProgress(), ' ', Percentage(), Bar()],
                          maxval=total).start()

    # the counter is only used for the optional progress bar
    for counter, pdb in enumerate(pdbs, 1):
        app.log.info("starting with PDB entry {0}.".format(pdb.upper()))

        # update the progress bar if required
        if args.progressbar: bar.update(counter)

        # test if the processed structure does already exist
        if args.incremental:
            
            # check if quaternary assemblies alreadt exist for this PDB entry
            if args.quat:
                
                # OEB files are divided into folders by 2nd & 3rd character of PDB code 
                assembly_dir = os.path.join(QUAT_OEB_DIR, pdb[1:3], pdb)
                
                # check if assembly directory contains OEB files 
                if len(glob.glob(os.path.join(assembly_dir, "{0}-*.oeb".format(pdb)))):
                    app.log.info("quaternary structure of {0} does already exist - skipped.".format(pdb))
                    continue
            
            # check if asymmetric unit structure exists already
            else:
                pass

        # path to the gzipped structure in the PDB mirror
        path = os.path.join(PDB_MIRROR_DIR, pdb[1:3], 'pdb{0}.ent.gz'.format(pdb))

        # check if pdb file exists on server
        if not os.path.exists(path):
            app.log.error("cannot read structure: PDB file for entry {0} does not exist.".format(path))
            continue

        ifs.open(str(path)) # must not be Unicode

        # by default, NMR models are treated as sequential molecules, which means
        # only the first model is considered here
        structure = ifs.GetOEGraphMols().next()
        structure.SetTitle(str(pdb)) # no Unicode here

        app.log.debug("structure loaded with {0} heavy atoms."
                      .format(OECount(structure, OEIsHeavy())))

        # extract information from the pdb header
        header = parse_header(structure)

        # debug PDB deposition date
        app.log.debug("structure deposition date is {0}.".format(header['REVDAT']))

        # ignore structures right now that are split
        if structure.GetBoolData('is_split'):
            app.log.warn("PDB Entry {0} is split over several PDB entries and "
                         "will be ignored.".format(pdb))
            continue

        OEAssignAromaticFlags(structure)

        # Determine hyb of all atoms in structure could be useful later on
        OEAssignHybridization(structure)

        # assign Bondi vdw radii
        OEAssignBondiVdWRadii(structure)

        # assign tripos atom names to all atoms
        OETriposAtomTypeNames(structure)

        # check for calpha structures and peptide ligands (short chains)
        CALPHA_FLAG = calc_calpha_ratio(structure)

        # ignore structure if at least one chain is below threshold
        if CALPHA_FLAG:
            app.log.warn("one of the polymer chains has a CA/CB ratio below the "
                         "minimum - structure will be skipped.")
            continue

        # set credo atom type flags to all atoms
        set_atom_type_flags(structure)
        
        # get information about ligands and polymers from the database
        pdb_polymer_info = db.get_pdb_polymer_info(pdb)
        pdb_ligand_info = db.get_pdb_ligand_info(pdb)
        pdb_sstruct_info = db.get_pdb_sstruct_info(pdb)

        # ASYMMETRIC UNIT DOES NOT CONTAIN ANY LIGANDS
        if not pdb_ligand_info:
            app.log.debug('assymmetric unit does not contain any ligands.')

        # assign entity serial numbers all atoms in ASU
        assign_entity_serials(structure, pdb_polymer_info, pdb_ligand_info,
                              pdb_sstruct_info)

        # identify the surface atoms of all entities in asymmetric structure
        identify_surface_atoms(structure)

        # generate quaternary assembly if PISA prediction is available
        if args.quat:

            biomolecules = []

            # check whether assemblies can be generated or if the asu has to split into monomers
            num_assemblies = db.get_pisa_num_assemblies(pdb)

            # no prediction in PISA - keep structure as it is
            if num_assemblies == -1:
                app.log.debug('no stable PISA prediction found.')
                biomolecules.append(structure)

            # structure is predicted to be monomeric
            elif num_assemblies == 0:
                app.log.debug('PISA predicts Monomer.')

                # check if structure has to be split into monomers, i.e. if it
                # has more than one chain
                pdb_chain_ids = set(OEAtomGetResidue(atom).GetChainID()
                                    for atom in structure.GetAtoms())

                # more than one chain / structure has to be split into monomers!
                if len(pdb_chain_ids) > 1:
                    app.log.debug("structure contains {0} chains but is predicted"
                                  " to be monomeric - structure will be split."
                                  .format(len(pdb_chain_ids)))

                    # iterate through chains and generate assemblies using the
                    # individual chains
                    for assembly_serial, pdb_chain_id in enumerate(sorted(pdb_chain_ids),1):
                        biomolecule = get_subsetmol(structure, OEHasChainID(pdb_chain_id))
                        biomolecule.SetIntData('assembly_serial', assembly_serial)

                        # get the ligand entity ids to attach only those that are
                        # actually found inside the new assembly
                        ligand_entity_serials = get_ligand_entity_serials(biomolecule)

                        # attach or remove from structure
                        if ligand_entity_serials:
                            biomolecule.SetData('ligand_entity_serials',
                                                ligand_entity_serials)
                        else:
                            biomolecule.DeleteData('ligand_entity_serials')

                        biomolecules.append(biomolecule)

                # STRUCTURE ONLY HAS A SINGLE CHAIN / NO SPLIT NECESSARY
                else: biomolecules.append(structure)

            # PISA contains predicted assemblies
            else:

                # get PISA data for quaternary assembly
                pisa = db.get_pisa_data(pdb)

                # generate quaternary assemblies using the PISA prediction
                biomolecules = biomol.generate_biomolecule(structure, pisa)
            
            # iterate through biomolecules and save structural data to disk 
            for biomolecule in biomolecules:

                # get the assembly number - important later on for the credo schema
                assembly_serial = biomolecule.GetIntData('assembly_serial')

                # bad PISA prediction
                if not biomolecule.NumAtoms():
                    app.log.warn("assembly {0} does not have any atoms!".format(assembly_serial))
                    continue

                # probably ligand-only structure (e.g. vancomycin)
                elif OECount(biomolecule, OEIsPolymerAtom) == 0:
                    app.log.info("assembly {0} does not have any polymer atoms."
                                 .format(assembly_serial))

                # get ligand molecules / ligand entity identifier are attached to
                # biomolecule
                ligands = get_ligands(biomolecule)
                
                # attach ligands as objects to biomolecule
                if ligands:
                    biomolecule.SetListData('ligands', ligands)
                    app.log.debug("{0} ligand(s) attached to the quaternary structure"
                                  " of assembly {1}.".format(len(ligands), assembly_serial))
                    
                    # save ligand structures on disk
                    # not really necessary because ligand molecules are attached
                    # to the parent molecule
                    if args.ligands: pass
                
                ### save structures to disk ###                
                
                # write assembly in openeye binary format
                if args.oeb:

                    # files are written to a divided file mirror
                    path = os.path.join(QUAT_OEB_DIR, pdb[1:3], pdb)
                    filename = os.path.join(path, '{0}-{1}.oeb'.format(pdb, assembly_serial))
                    
                    # make new directories recursively if necessary
                    if not os.path.exists(path): os.makedirs(path)                                   
                    
                    # write the molecule
                    oebofs.open(str(filename)) # no Unicode!
                    OEWriteMolecule(oebofs, biomolecule)

                    app.log.debug("biological assembly {0} written as OEB to {1}"
                                  .format(assembly_serial, filename))
                
                # write biomolecule as PDB file for easier visual inspection -
                # warning: can produce non-PDB format compliant files
                if args.pdb:
                    
                    # files are written to a divided file mirror
                    path = os.path.join(QUAT_PDB_DIR, pdb[1:3], pdb)
                    filename = os.path.join(path, '{0}-{1}.pdb'.format(pdb, assembly_serial))
                    
                    # make new directories recursively if necessary
                    if not os.path.exists(path): os.makedirs(path)                                   
                    
                    # write the molecule
                    pdbofs.open(str(filename)) # no Unicode!
                    OEWriteMolecule(pdbofs, biomolecule)

                    app.log.debug('biological assembly {0} written as PDB to {1}'
                                  .format(assembly_serial, filename))
                
                # write chain structures to disk
                if args.chains:

                    # get all chain ids
                    pdb_chain_ids = set(OEAtomGetResidue(atom).GetChainID()
                                        for atom in biomolecule.GetAtoms())

                    # iterate through all chains of biomolecule
                    for pdb_chain_id in pdb_chain_ids:

                        # only keep polymer sequence / ignore ligands and solvents
                        predicate = OEAndAtom(OEHasChainID(pdb_chain_id),
                                              OEIsPolymerAtom)

                        # create subset molecule for the complete chain
                        chain = get_subsetmol(biomolecule, predicate)

                        # chain probably contains only non-polymer residues
                        if chain.NumAtoms() == 0:
                            app.log.debug("skipped saving of chain {0} of biomolecule"
                                          " {1} because it did not contain any polymers."
                                          .format(pdb_chain_id, assembly_serial))
                            continue

                        # remove PDB header
                        OEClearPDBData(chain)

                        # remove ligands from chain
                        chain.DeleteData('ligands')

                        # set chain identifying data
                        chain.SetStringData('pdb', pdbcode)
                        biomolecule.SetIntData('assembly_serial', assembly_serial)
                        chain.SetStringData('pdb_chain_id', pdb_chain_id)

                        # write chain to disk
                        filename = '{0}_{1}_{2}.oeb'.format(pdbcode, assembly_serial,
                                                            pdb_chain_id)
                        
                        path = os.path.join(QUAT_OEB_CHAIN_DIR, pdb[1:3], pdb,
                                            filename)

                        # make new directories recursively if necessary
                        if not os.path.exists(path): os.makedirs(path)

                        # opene binary output stream
                        oebofs.open(str(path))
                        OEWriteMolecule(oebofs, biomolecule)

                        # debug chain saving info
                        app.log.debug("chain {0} of biomolecule {1} with {2} atoms"
                                      " written to {3}.".format(pdb_chain_id,
                                                                assembly_serial,
                                                                chain.NumAtoms(),
                                                                path))
        
        app.log.debug("finished processing PDB entry {0}.".format(pdb.upper()))
    
    # finish the progress bar
    if args.progressbar: bar.finish()  