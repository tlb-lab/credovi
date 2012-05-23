"""
python credovi.py mmcif currentpdbs | parallel --eta --halt 2 -n 4 python credovi.py credo contacts -Q -S {1},{2},{3},{4}
"""

import os
import csv
import glob
import string
from math import sqrt # faster than the numpy version

from eyesopen.oechem import *

from credovi import app
from credovi.structbio import structure as struct
from credovi.structbio.interactions import (is_hbond, is_weak_hbond, is_xbond,
                                            is_aromatic, is_carbonyl, is_metal_complex,
                                            is_ionic, is_hydrophobic)

from credovi.util.timer import Timer

# register new dialect to write tab-delimited files
csv.register_dialect('tabs', delimiter='\t')

# directory containing all the biological assemblies in OEB format
QUAT_OEB_DIR = app.config.get('directories','quat_oeb')

def get_assembly_sets(args):
    """
    Returns a list of assembly sets in the form [(PDB, [1XXX-0.oeb, 1XXX-1.oeb,...]),...].
    Only used for cryst. All assemblies must have been processed before.
    """
    def _get_assembly_sets(pdbs):
        """
        """
        assemblysets = []

        # iterate through all specified pdb codes / remove possible newline character
        for pdb in pdbs:

            # keep all assemblies of the set together
            assemblyset = []

            # OEB files are divided into folders by 2nd & 3rd character of PDB code
            assembly_dir = os.path.join(QUAT_OEB_DIR, pdb[1:3].lower(), pdb.lower())

            if os.path.exists(assembly_dir):

                # add all assemblies in directory to our assembly set
                for oeb in os.listdir(assembly_dir): assemblyset.append(oeb)

                # add assembly set to all the other assembly sets
                if assemblyset: assemblysets.append((pdb.upper(), assemblyset))

            # no assemblies exists for the given PDB code
            else:
                app.log.warn("no assembly set found for PDB entry: {pdb}".format(pdb=pdb))

        return assemblysets

    # set mirror to quaternary assemblies
    if args.quat:

        # get all available assemblies
        if args.all:
            pdbs = []

        # a list of PDB codes was provided on the command line
        elif args.structures:
            pdbs = [pdb.upper() for pdb in args.structures.strip().split(',') if pdb]

        # PDB codes are taken from a test set in the config file
        elif args.testset:
            pdbs = map(string.upper, app.config['test sets'][args.testset])

        assemblysets = _get_assembly_sets(pdbs)

    # get only the asymmetric unit structures
    else:
        assemblysets = []

    return sorted(assemblysets)

def write_ligands(ligands, pdb, biomolecule, writer):
    """
    """
    for ligand in ligands:
        res_num = ligand.GetIntData('res_num') if ligand.HasData('res_num') else '\N'

        row = [pdb, biomolecule, ligand.GetIntData('entity_serial'),
               ligand.GetStringData('pdb_chain_id'), res_num,
               ligand.GetStringData('name'), OECount(ligand, OEIsHeavy())]

        writer.writerow(row)

def write_aromatic_rings(aromatic_rings, pdb, biomolecule, writer):
    """
    """
    aromatic_ring_serial = 0
    for aromatic_ring_serial, atoms in aromatic_rings:
        for atom in atoms:
            atom_serial = OEAtomGetResidue(atom).GetSerialNumber()

            row = [pdb, biomolecule, aromatic_ring_serial, atom_serial]

            writer.writerow(row)

    app.log.debug("{0} aromatic ring systems found.".format(aromatic_ring_serial))

def write_residues(structure, biomolecule, residues, writer):
    """
    """
    for residue in residues:
        for atom in structure.GetAtoms(OEAtomIsInResidue(residue)):

            # residue data is atom specific therefore the residue has to be obtained
            # for every residue atom
            atomres = OEAtomGetResidue(atom)

            row = ['\N' for i in range(32)]

            # PDB DATA
            row[0]  = structure.GetTitle() # PDB
            row[1]  = biomolecule
            row[2]  = atomres.GetSerialNumber()
            row[3]  = 'HETATM' if atomres.IsHetAtom() else 'ATOM'
            row[4]  = atom.GetName().strip()
            row[5]  = atomres.GetAlternateLocation()
            row[6]  = atomres.GetName().strip()
            row[7]  = atomres.GetChainID()
            row[8]  = atom.GetStringData('pdb_chain_asu_id') if atom.GetStringData('pdb_chain_asu_id') else atomres.GetChainID()
            row[9]  = atomres.GetResidueNumber()
            row[10] = atomres.GetInsertCode()
            row[11] = '{%s,%s,%s}' % structure.GetCoords(atom) # format required for vector3d extension
            row[12] = atomres.GetOccupancy()
            row[13] = atomres.GetBFactor()
            row[14] = OEGetAtomicSymbol(atom.GetAtomicNum()) and OEGetAtomicSymbol(atom.GetAtomicNum()) or '\N'

            # ATOM HYBRIDISATION
            row[15] = atom.GetHyb()

            # TRIPOS ATOM TYPE
            row[16] = atom.GetType()

            # ENTITY_ID
            row[17] = atom.GetIntData('entity_serial')

            # ENTITY TYPE BITMASK
            row[18] = atom.GetIntData('entity_type_bm')

            # ATOM FLAGS
            row[19] = atom.GetIntData('hbond donor')            # DONOR
            row[20] = atom.GetIntData('hbond acceptor')         # ACCEPTOR
            row[21] = 1 if atom.IsAromatic() else 0             # AROMATIC
            row[22] = atom.GetIntData('weak hbond acceptor')    # WEAK_ACCEPTOR
            row[23] = atom.GetIntData('weak hbond donor')       # WEAK_DONOR
            row[24] = atom.GetIntData('hydrophobe')             # HYDROPHOBE
            row[25] = 1 if atom.IsMetal() else 0                # METAL
            row[26] = atom.GetIntData('pos ionisable')          # POS_IONIZABLE
            row[27] = atom.GetIntData('neg ionisable')          # NEG_IONIZABLE
            row[28] = atom.GetIntData('xbond donor')            # HALOGEN_BOND_DONOR
            row[29] = atom.GetIntData('xbond acceptor')         # HALOGEN_BOND_ACCEPTOR
            row[30] = atom.GetIntData('carbonyl oxygen')        # CARBONYL OXYGEN
            row[31] = atom.GetIntData('carbonyl carbon')        # CARBONYL CARBON

            writer.writerow(row)

def do(controller):
    """
    """
    # get the controller command
    cmd = controller.command

    # get the command line arguments and options
    args = controller.pargs

    # timer to clock functions and parts of the program
    timer = Timer()
    timer.start("app")

    # get the files that are to be processed
    assemblysets = get_assembly_sets(args)

    # iterate through all assembly sets that are going to be processed
    for pdb, assemblyset in assemblysets:
        app.log.info("starting with PDB entry {0}.".format(pdb))
        timer.start('assemblyset')

        # create a path for this structure to which all data will be written
        # PDB 'divided' directory convention is used here
        struct_data_dir = os.path.join(app.config.get('directories','credo_data'),
                                       pdb[1:3].lower(), pdb.lower())

        # make necessary directories recursively if not existing yet
        if not os.path.exists(struct_data_dir):
            os.makedirs(struct_data_dir)

            app.log.debug("created new directory {0} for PDB entry {1}."
                          .format(struct_data_dir, pdb))

        # open file handles for each PDB entry
        # too many filehandles for with statement...
        atomfs = open(os.path.join(struct_data_dir, 'atoms.credo'), 'w')
        contactfs = open(os.path.join(struct_data_dir, 'contacts.credo'), 'w')
        ligandfs = open(os.path.join(struct_data_dir, 'ligands.credo'), 'w')
        aromaticringfs = open(os.path.join(struct_data_dir, 'aromaticrings.credo'), 'w')
        chainfs = open(os.path.join(struct_data_dir, 'chains.credo'), 'w')

        # CSV writers
        ligwriter = csv.writer(ligandfs, dialect='tabs')
        aromaticringwriter = csv.writer(aromaticringfs, dialect='tabs')
        atomwriter = csv.writer(atomfs, dialect='tabs')
        cswriter = csv.writer(contactfs, dialect='tabs')
        chainwriter = csv.writer(chainfs, dialect='tabs')

        # iterate through each assembly of set and generate credo data for each
        for assembly in assemblyset:

            # we are only interested in the filename here
            filename, extension = os.path.splitext(assembly)
            timer.start('assembly')

            # biological assemblies already exist and are simply loaded from the directory
            if args.quat:
                pdbcode, biomolecule = filename.split('-')
                path = os.path.join(QUAT_OEB_DIR, pdb[1:3].lower(), pdb.lower(), assembly)

            # only asymmetric unit structures are supposed to be processed
            else: pass

            # this should not happen if assemblysets were loaded with get_assembly_sets()
            if not os.path.isfile(path):
                app.log.error("cannot process {0}: file {1} does not exist!".format(pdb, path))
                continue

            # load the structure as OEGraphMol
            # returns None if structure could not be loaded
            structure = struct.get_structure(str(path), str(pdb))

            # unable to parse structure: must be dodgy file
            if not structure:
                app.log.error("cannot process {0}: unable to parse {1}!".format(pdb, path))
                continue

            # debug number of atoms in structure
            app.log.debug("loaded structure with {0} atoms.".format(structure.NumAtoms()))

            # Hybridization attribute does not seem to be preserved in OEB
            OEAssignHybridization(structure)

            # assembly serial number
            assembly_serial = structure.GetIntData('assembly_serial')

            app.log.debug("structure is assembly number {0}.".format(assembly_serial))

            # get the ligands from the structure / ligands are already attached
            ligands = structure.GetListData('ligands')

            app.log.debug("{0} ligand(s) found in structure.".format(len(ligands)))

            # write ligand information to separate file
            write_ligands(ligands, pdb, biomolecule, ligwriter)

            # get all aromatic rings in the structure
            aromatic_rings = struct.get_aromatic_rings(structure,
                                                       app.config['atom types']['aromatic'].values())

            # write aromatic rings to file
            write_aromatic_rings(aromatic_rings, pdb, biomolecule, aromaticringwriter)

            timer.start()

            ### get all by all contacts / includes covalently bound neighbour atoms
            contacts = OEGetNearestNbrs(structure, app.config['cutoffs']['cutoff'],
                                        OENearestNbrsMethod_Auto)

            # debug how much time it took to get all contacts
            app.log.debug("all contacts identified in {0:.2f} seconds.".format(timer.elapsed()))

            # protonate structure if necessary in order to calculate hbond angles
            if not OEHasExplicitHydrogens(structure):
                OEAddExplicitHydrogens(structure, False, True)
                OESet3DHydrogenGeom(structure)

            timer.start()

            # keep track of all residues that are interacting and need to be in credo
            interacting_residues = set()

            for contact in contacts:
                """
                SIFT:
                    0: CLASH
                    1: COVALENT
                    2: VDW_CLASH
                    3: VDW
                    4: PROXIMAL
                    5: HBOND
                    6: WEAK_HBOND
                    7: HALOGEN_BOND
                    8: IONIC
                    9: METAL_COMPLEX
                   10: AROMATIC
                   11: HYDROPHOBIC
                   12: CARBONYL
                """
                # initialize structural interaction fingerprint
                SIFt = [0] * 13

                IS_INTRAMOLECULAR = False

                # get interacting atoms
                if contact.GetBgn() < contact.GetEnd(): atom_bgn, atom_end = contact.GetBgn(), contact.GetEnd()
                else: atom_bgn, atom_end = contact.GetEnd(), contact.GetBgn()

                # first atom is entity atom
                if atom_bgn.GetIntData('entity_serial') > 0:

                    # ignore contacts of non-exposed entity atoms
                    if atom_bgn.GetIntData('is_exposed') == 0: continue

                    # second atom is entity atom
                    if atom_end.GetIntData('entity_serial') > 0:

                        # ignore contacts of non-exposed entity atoms
                        if atom_end.GetIntData('is_exposed') == 0: continue

                        # set a flag for intra-entity contacts
                        if atom_bgn.GetIntData('entity_serial') == atom_end.GetIntData('entity_serial'):
                            IS_INTRAMOLECULAR = True

                # first atom is solvent
                else:

                    # second atom is solvent / ignore solvent-solvent contacts
                    if atom_end.GetIntData('entity_serial') == 0: continue

                    # second atom belongs to entity / ignore contacts with non-exposed entity atoms
                    if atom_end.GetIntData('is_exposed') == 0: continue

                # get the residues
                res_bgn, res_end = OEAtomGetResidue(atom_bgn), OEAtomGetResidue(atom_end)

                # ignore all intra-residue contacts
                if OESameResidue(res_bgn, res_end): continue

                # get interatomic distance
                distance = sqrt(contact.GetDist2())

                # ignore intra-molecular contacts that are only separated by three bonds
                if IS_INTRAMOLECULAR and OEGetPathLength(atom_bgn, atom_end, 2):
                    continue

                # get the sum of van der waals radii of both atoms
                sum_vdw_radii = atom_bgn.GetRadius() + atom_end.GetRadius()

                # get the sum of covalent radii
                sum_cov_radii = OEGetCovalentRadius(atom_bgn.GetAtomicNum()) + OEGetCovalentRadius(atom_end.GetAtomicNum())

                # use the connection table to identify covalent bonds
                # can identify covalent bonds of unusual length as well
                if atom_end in atom_bgn.GetAtoms(): SIFt[1] = 1

                # check for atomic clash of non-bonded atoms
                elif distance < sum_cov_radii - 0.1: SIFt[0] = 1

                # check vdw radii
                elif distance < sum_vdw_radii: SIFt[2] = 1
                elif distance <= sum_vdw_radii + 0.1: SIFt[3] = 1 # + vdw comp factor

                # label as proximal
                else: SIFt[4] = 1

                # skip this step if atoms are clashing, covalently bonded or above
                # the maximum feature contact type distance (4.5) to save time
                if not any(SIFt[:2]) and distance < app.config['cutoffs']['contact type dist max']:

                    SIFt[5]  = is_hbond(structure, atom_bgn, atom_end, distance)
                    SIFt[6]  = is_weak_hbond(structure, atom_bgn, atom_end, distance)
                    SIFt[7]  = is_xbond(structure, atom_bgn, atom_end, distance, sum_vdw_radii)
                    SIFt[8]  = is_ionic(atom_bgn, atom_end, distance)
                    SIFt[9]  = is_metal_complex(atom_bgn, atom_end, distance)
                    SIFt[10] = is_aromatic(atom_bgn, atom_end, distance)
                    SIFt[11] = is_hydrophobic(atom_bgn, atom_end, distance)
                    SIFt[12] = is_carbonyl(atom_bgn, atom_end, distance)

                ### write contact details to file / atom serial is used to identify atom

                row = ['\N' for i in range(7)]

                row[0] = pdb # PDB
                row[1] = biomolecule
                row[2] = res_bgn.GetSerialNumber()
                row[3] = res_end.GetSerialNumber()
                row[4] = distance

                # get the structural interaction type for this inter-atomic contact
                # as the combination of both individual entity type bit masks /
                # the first bit mask will be shifted by 6 positions to make space
                # for the second / this way directionality is kept and the integer
                # size below 4096
                row[5] = (atom_bgn.GetIntData('entity_type_bm') << 6) + atom_end.GetIntData('entity_type_bm')

                # write boolean as 0 or 1
                row[6] = 1 if IS_INTRAMOLECULAR else 0

                row.extend(SIFt)

                cswriter.writerow(row)

                # add both residues to the set of interacting residues
                # their atoms will be written to a file later on
                interacting_residues.add(res_bgn)
                interacting_residues.add(res_end)

            app.log.debug("all contacts processed in {0:.2f} seconds."
                          .format(timer.elapsed()))

            # remove explicit hydrogens again from structure before writing atoms
            OESuppressHydrogens(structure)

            timer.start()

            # write all residue atoms interacting with different entities to file
            write_residues(structure, biomolecule, interacting_residues, atomwriter)

            # time it took to write all residue atoms to a filestream
            app.log.debug("all atoms written in {0:.2f} seconds.".format(timer.elapsed()))

            # time passed to finish processing the complete assembly
            app.log.debug("finished processing assembly in {0:.2f} seconds."
                          .format(timer.elapsed('assembly')))

            atomfs.flush()
            contactfs.flush()
            ligandfs.flush()
            aromaticringfs.flush()
            chainfs.flush()

        # close files after assembly set has been processed
        atomfs.close()
        contactfs.close()
        ligandfs.close()
        aromaticringfs.close()
        chainfs.close()

        app.log.info("finished processing assembly set for PDB entry {0} in "
                     "{1:.2f} seconds.".format(pdb, timer.elapsed('assemblyset')))

    app.log.info("finished processing all structures in {0} days/h/m/s."
                 .format(timer.formatted('app')))