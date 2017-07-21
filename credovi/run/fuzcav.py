from __future__ import absolute_import

import numpy as np
from sqlalchemy import and_, func, or_, Table
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

from credoscript import Session, metadata
from credoscript.models import Atom, Ligand, Peptide, BindingSiteResidue

from credovi import app
from credovi.schema import engine
from credovi.util.timer import Timer
from credovi.schema.tables.ligands import binding_site_fuzcav

from credovi.structbio import fuzcav

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

    insert = binding_site_fuzcav.insert()
    tracker = fuzcav.get_tracker()

    # get the fuzcav side chain representative table from the credoscript metadata
    metadata.reflect(schema='bio', only=('fuzcav_rep_sc_atoms',))
    fuzcav_rep_sc_atoms = Table('bio.fuzcav_rep_sc_atoms', metadata, autoload=True)

    timer.start()

    session = Session()

    # get all ligands that have more than 7 heavy atoms and no clashes
    query = session.query(Ligand.ligand_id, Ligand.biomolecule_id)
    query = query.filter(and_(Ligand.num_hvy_atoms>=7, Ligand.is_clashing==False))

    if args.incremental:

        # subquery to get the current max ligand_id from the binding_site_fuzcav table
        sq = session.query(func.max(binding_site_fuzcav.c.ligand_id).label('ligand_id')).subquery('sq')

        # only include new ligands
        query = query.filter(Ligand.ligand_id > sq.c.ligand_id)

    ligand_ids = query.order_by(Ligand.ligand_id).all()

    # debug how much time it took to get all contacts
    app.log.debug("all new ligand identifiers retrieved in {0:.2f} seconds."
                  .format(timer.elapsed()))

    #
    query = BindingSiteResidue.query.join('Peptide', 'Atoms')
    #query = query.join(Peptide, Peptide.residue_id==BindingSiteResidue.residue_id)
    #query = query.join(Atom, Atom.residue_id==Peptide.residue_id)
    query = query.outerjoin(fuzcav_rep_sc_atoms, and_(fuzcav_rep_sc_atoms.c.res_name==Peptide.res_name,
                                                      fuzcav_rep_sc_atoms.c.atom_name==Atom.atom_name))
    query = query.filter(and_(Peptide.is_non_std==False,
                              or_(Atom.atom_name=='CA',
                                  fuzcav_rep_sc_atoms.c.atom_name!=None)))
    query = query.with_entities(Peptide.res_name, Atom)

    if args.progressbar:
        bar = ProgressBar(widgets=['Binding Sites: ', SimpleProgress(), ' ',
                                   Percentage(), Bar()], maxval=len(ligand_ids)).start()

    # iterate through ligands
    for counter, row in enumerate(ligand_ids, 1):
        if args.progressbar: bar.update(counter)
        ligand_id, biomolecule_id = row.ligand_id, row.biomolecule_id

        timer.start()

        # get all the fuzcav atoms (either CA or representative)
        # important to use the proper atom partition!
        atoms = query.filter(and_(BindingSiteResidue.ligand_id==ligand_id,
                                  Atom.biomolecule_id==biomolecule_id)).all()

        # debug how much time it took to get all contacts
        app.log.debug("all FuzCav atoms retrieved in {0:.2f} seconds."
                      .format(timer.elapsed()))

        # ignore hits with too few peptides
        if len(atoms) < 14:
            app.log.debug("Ligand {} has only {} FuzCav atoms and will be "
                          "ignored.".format(ligand_id, len(atoms)))
            continue

        # get the calpha atom and its features for each residue
        calphas = ((np.array(atom.coords, dtype=float), (fuzcav.FEATURES[res_name]))
                   for res_name, atom in atoms if atom.atom_name=='CA')

        # get the representative atom and its features for each residue
        representatives = ((np.array(atom.coords, dtype=float), (fuzcav.FEATURES[res_name]))
                           for res_name, atom in atoms
                           if atom.atom_name==fuzcav.REPRESENTATIVES[res_name])

        timer.start()

        calphafp = fuzcav.make_fp(calphas, tracker)
        repfp = fuzcav.make_fp(representatives, tracker)

        # debug how much time it took to get all contacts
        app.log.debug("fingerprints generated in {0:.2f} seconds."
                      .format(timer.elapsed()))

        # insert the fingerprints into the table
        if not args.dry_run:
            engine.execute(insert, ligand_id=ligand_id,
                           calphafp=calphafp.tolist(), repfp=repfp.tolist())

    # finish the optional progress bar
    if args.progressbar: bar.finish()

    session.close()
