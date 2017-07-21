from itertools import groupby, ifilter
from operator import itemgetter
from string import ascii_uppercase

from eyesopen.oechem import *

from credovi import app
from credovi.util.timer import timer
from credovi.structbio.structure import get_ligand_entity_serials

# dictionary of possible chain ids and priority number
CHAIN_IDS       = {'1':53, '0':52, '3':55, '2':54, '5':57, '4':56, '7':59, '6':58,
                   '9':61, '8':60, 'A':0, 'C':2, 'B':1, 'E':4, 'D':3, 'G':6, 'F':5,
                   'I':8, 'H':7, 'K':10, 'J':9, 'M':12, 'L':11, 'O':14, 'N':13,
                   'Q':16, 'P':15, 'S':18, 'R':17, 'U':20, 'T':19, 'W':22, 'V':21,
                   'Y':24, 'X':23, 'Z':25, 'a':26, 'c':28, 'b':27, 'e':30, 'd':29,
                   'g':32, 'f':31, 'i':34, 'h':33, 'k':36, 'j':35, 'm':38, 'l':37,
                   'o':40, 'n':39, 'q':42, 'p':41, 's':44, 'r':43, 'u':46, 't':45,
                   'w':48, 'v':47, 'y':50, 'x':49, 'z':51}

@timer("water clashes removed in {0:.2f} seconds.")
def remove_water_clashes(structure):
    """
    Removes all water atoms from the structure that are clashing (happens usually
    after biomolecule transformations).
    """
    # Identify clashing water atoms // o-h...o distance is 1.97
    clashes = OEGetNearestNbrs(structure, 2.0, OENearestNbrsMethod_Auto)

    # keep only clashing water atoms
    clashes = set([clash.GetEnd() for clash in clashes if
                   OEGetResidueIndex(clash.GetBgn())==OEResidueIndex_HOH
                   and OEGetResidueIndex(clash.GetEnd())==OEResidueIndex_HOH])

    # delete clashing atoms from structure
    for water in clashes: structure.DeleteAtom(water)

@timer("biological assemblies were generated in {0:.2f} seconds.")
def generate_biomolecule(structure, biomt):
    """
    Generates the quaternary assembly of a PDB structure using the top-ranked
    PISA prediction (if available).
    """
    # get the PDB code here to prevent concatenated codes as titles
    pdbcode = structure.GetTitle()

    biomolecules = []

    # an assembly equates to a valid biomolecule
    # a PDB structure can contain more than one biomolecule in the asu (1c3h)
    for assembly_serial in sorted(biomt):
        biomolecule = OEMol()

        # preserve original PDB header
        OECopyPDBData(biomolecule, structure)

        # dictionary to keep track of the transformations that were executed
        # on newly created chains
        biomol_trans_map = {}

        #
        cur_chain_id = None

        # list of all possible new pdb chain identifiers ordered by priority (
        # A-Z,a-z,0-9)
        new_chain_ids = [pdb_chain_id for pdb_chain_id, priority in
                         sorted(CHAIN_IDS.items(), key=itemgetter(1))]

        # remove existing chain identifiers - will be used for identity operations
        for pdb_chain_id in set(OEAtomGetResidue(atom).GetChainID()
                                for atom in structure.GetAtoms()):
            new_chain_ids.remove(pdb_chain_id)

        # get the largest entity serial number in the ASU
        # new entities generated through symmetry operations have to above that
        struct_max_entity_serial = max(atom.GetIntData('entity_serial')
                                       for atom in structure.GetAtoms())

        # iterate through chains that are going to be transformed
        for pdb_chain_id in biomt[assembly_serial]:

            # chain of original structure that is going to be transformed
            try:
                molecule = get_subsetmol(structure,OEHasChainID(pdb_chain_id))

            # The PDB chain ID is too long - problably parsing error of REMARK350
            except TypeError:
                app.log.warn("BIOMT PDB chain identifier {} is not a single character!"
                             .format(pdb_chain_id))
                continue

            # somehow PISA has chains that are not in the pdb entry:
            if not molecule.NumAtoms():
                app.log.warn("BIOMT PDB chain identifier {0} cannot be found "
                             "in structure!".format(pdb_chain_id))
                continue

            # remove PDB header from chain
            OEClearPDBData(molecule)

            # iterate through biounit transformations
            for operation_serial, details in biomt[assembly_serial][pdb_chain_id].items():

                # copy of the object that will be transformed
                entity = OEMol(molecule)

                # transformation details
                rotation = details['rotation']
                translation = details['translation']
                is_at_identity = details['is_at_identity']

                ### transform entity and add to assembly ###

                # chain is at identity, no transformation necessary
                if is_at_identity:
                    msg = "identity matrix found for assembly {0} operation {1} chain {2}."
                    msg = msg.format(assembly_serial, operation_serial, pdb_chain_id)

                # molecule has to be transformed
                else:

                    # transformation can only be done on conformer
                    conformer = entity.GetConf(OEHasConfIdx(0))

                    # ROTATE AND TRANSLATE (IN THIS ORDER!)
                    OERotMatrix(rotation).Transform(conformer)
                    OETranslation(translation).Transform(conformer)

                    msg = "Chain {0} rotated by {1} and translated by {2} "
                    msg = msg.format(pdb_chain_id,
                                     "%4.2f " * len(rotation) % tuple(rotation),
                                     "%5.2f " * len(translation) % tuple(translation))

                # first operation on assembly
                if operation_serial == 1: app.log.debug(msg)

                # not the first operation - new chain id and entity ids are required
                else:

                    # get new chain identifier, try to keep original
                    cur_chain_id = new_chain_ids.pop(0)

                    # keep track of the transformations that were performed
                    # on this chain
                    biomol_trans_map[cur_chain_id] = details

                    # debug information about the new pdb chain id that was created
                    app.log.debug(msg + "to generate chain {0}.".format(cur_chain_id))

                    # set chain and entity id for all atoms
                    for atom in entity.GetAtoms():

                        # set new chain id
                        residue = OEAtomGetResidue(atom)
                        residue.SetChainID(cur_chain_id)

                        # label atom as quaternary and keep track of old pdb chain id
                        # important later on to map external data
                        atom.SetIntData('is_quaternary',1)
                        atom.SetStringData('pdb_chain_asu_id', pdb_chain_id)

                        # only increase entity serial number for non-solvents
                        if OENotIsWater(atom):

                            # increment entity serial number
                            new_entity_serial = atom.GetIntData('entity_serial') + struct_max_entity_serial * operation_serial
                            atom.SetIntData('entity_serial', new_entity_serial)

                #


                # add newly created chain to biomolecule
                OEAddMols(biomolecule, entity)

        #

        # remove water molecules that clash in quaternary assembly
        remove_water_clashes(biomolecule)

        # get all the ligand entity serial numbers from the current assembly
        ligand_entity_serials = get_ligand_entity_serials(biomolecule)

        # attach or remove from assembly
        if ligand_entity_serials:
            biomolecule.SetData('ligand_entity_serials', ligand_entity_serials)
        else:
            biomolecule.DeleteData('ligand_entity_serials')

        # debug ligand entity serial number info in biomolecule
        app.log.debug("Ligand entity serial numbers in assembly {0}: {1} ({2})"
                      .format(assembly_serial, ligand_entity_serials, biomolecule.HasData('ligand_entity_serials')))

        biomolecule.SetTitle(pdbcode)
        biomolecule.SetIntData('assembly_serial', assembly_serial)

        # identify alternate location groups in biomolecule and assign new alternate location codes
        # or remove completely superimposed residues
        set_alternate_locations(biomolecule)

        # reorder the atoms using the residue information associated with each atom
        OEPDBOrderAtoms(biomolecule)

        # assign PDB serial numbers to new atoms
        OEAssignSerialNumbers(biomolecule)

        #
        OEPerceiveSecondaryStructure(biomolecule)

        biomolecules.append(biomolecule)

    return biomolecules

@timer("alternate locations set for disordered ligands in {0:.2f} seconds.")
def set_alternate_locations(structure):
    '''
    This function looks for residues that might be disordered after biomolecule
    transformations were applied. It identifies clashing residues and calculates
    the RMSD to decide whether the residues in close contact are superimposed and
    one of them can be discarded or if a alternation location group has to be created.
    The function also merges the entity_serial and pdb_chain_asu_id information
    of the two residues.
    '''
    # only consider ligands for alternate location group analysis
    if not structure.HasData('ligand_entity_serials'):
        app.log.debug("Alternate location group analysis for biomolecule {0} "
                      "skipped because no ligands were found."
                      .format(structure.GetIntData('assembly_serial')))

        return

    # some entity serial numbers become obsolete if they are superimposed after
    # transformation, so keep track of them right now
    entity_serials = set(structure.GetData('ligand_entity_serials'))

    # get residues that are likely to be disordered in the crystal
    clashes = OEGetNearestNbrs(structure, 2, OENearestNbrsMethod_Auto)
    clashes = ((cs.GetBgn(), cs.GetEnd()) for cs in clashes)

    # only include clashes between ligand residues
    clashes = ifilter(lambda (bgn,end): (bgn.GetIntData('entity_type_bm') == 2 and end.GetIntData('entity_type_bm') == 2), clashes)

    # also only consider newly created quaternary atoms
    clashes = ifilter(lambda (bgn,end): (bgn.GetIntData('is_quaternary') or end.GetIntData('is_quaternary')), clashes)

    # only keep clashes between the same atoms (equal pdb names)
    clashes = ifilter(lambda (bgn,end): bgn.GetName() == end.GetName(), clashes)

    # get possibly disordered residues
    clashes = set((OEAtomGetResidue(bgn), OEAtomGetResidue(end)) for bgn,end in clashes)

    # keep only residues with the same name
    clashes = ifilter(lambda (x,y): x.GetName()==y.GetName() and not OESameResidue(x,y), clashes)

    # and residue number (1C6L)
    clashes = ifilter(lambda (x,y): x.GetResidueNumber()==y.GetResidueNumber(), clashes)

    # create alternate location groups
    map = {}
    for res_bgn, res_end in clashes:
        if res_bgn in map: map[res_bgn].add(res_end)
        else: map[res_bgn] = set([res_bgn, res_end])

        if res_end in map: map[res_end].add(res_bgn)
        else: map[res_end] = set([res_bgn,res_end])

    # compile a list of possible alternate groups
    # residues are sorted by pdb identifier so asymmetric pdb chain ids should come first
    altgroups = set(tuple(sorted(row)) for row in map.values())

    app.log.debug("{0} possible alternate location groups found in assembly."
                  .format(len(altgroups)))

    # proceed to either deleting superimposed ligands or setting new alternate locations
    if altgroups:

        # calculate the RMSD between residues in altgroups to decide whether
        # residue can be removed / use list() to avoid error if the set size is
        # changed during iteration
        for altgroup in list(altgroups):
            refres = altgroup[0]
            refmol = get_subsetmol(structure, OEIsResidue(refres))

            # fitres should always be quat because altgroup is sorted by pdb chain id
            for fitres in altgroup[1:]:
                fitmol = get_subsetmol(structure, OEIsResidue(fitres))

                # calculate RMSD
                rmsd = OERMSD(refmol, fitmol, True, True, False)

                # delete superimposed residue from structure
                if rmsd < app.config['cutoffs']['RMSD disorder cutoff']:

                    # get all residue atoms
                    for atom in structure.GetAtoms(OEAtomIsInResidue(fitres)):

                        # remove atom as well as entity id
                        entity_serials.discard(atom.GetIntData('entity_serial'))
                        structure.DeleteAtom(atom)

                    # discard alternate group because superimposed residues were
                    # deleted
                    altgroups.discard(altgroup)

                    app.log.debug("quaternary {0} with RMSD {1:.2f}A removed."
                                  .format(fitres,rmsd))

                else:

                    # throw a warning if the residues have a high rmsd
                    # most likely only very atoms are actually clashing and where
                    # disorderd before rotation
                    if rmsd > 3.0:
                        app.log.warn("quaternary residue {0} found with RMSD {1:.2f}A"
                                     .format(fitres,rmsd))

        ### merge residue information of alternate location groups by taking the
        ### first residue as a reference to make sure that entity serial number
        ### and pdb chain id are consistent

        # iterate through all identified alternation location groups in this biomolecule
        for altgroup in altgroups:

            # list of available alternate location codes
            avl_altlocs = [c for c in ascii_uppercase]

            # use the first residue as reference
            residue_ref = altgroup[0]

            # get the first pdb chain id of altgroup as reference
            pdb_chain_id_ref = residue_ref.GetChainID()

            # get the entity id and pdb chain asu information from the first atom
            refatom = structure.GetAtoms(OEAtomIsInResidue(residue_ref)).next()
            entity_id_ref = refatom.GetIntData('entity_serial')
            pdb_chain_asu_id_ref = refatom.GetStringData('pdb_chain_asu_id')

            # check if altgroup already contains alternate locations and split further
            # if necessary (see 1DVX B 125 DIF)
            for residue in altgroup:

                # sort residue atoms by alternate location
                atoms = sorted(structure.GetAtoms(OEAtomIsInResidue(residue)), key=lambda x: OEAtomGetResidue(x).GetAlternateLocation())

                # group residue atom by alternate location code in order to
                # increment them if necessary
                for altloc, atomiter in groupby(atoms, key=lambda x: OEAtomGetResidue(x).GetAlternateLocation()):

                    # create new unique alternate location for this group
                    altloc_new = avl_altlocs.pop(0)

                    for atom in atomiter:

                        atomres = OEAtomGetResidue(atom)

                        # set the reference chain id so every altloc has the same
                        atomres.SetChainID(pdb_chain_id_ref)

                        # set the new alternate location
                        atomres.SetAlternateLocation(altloc_new)

                        # remove obsolete entity id
                        entity_id_altgroup = atom.GetIntData('entity_serial')

                        # only remove is residue is not the reference
                        if entity_id_altgroup != entity_id_ref: entity_serials.discard(entity_id_altgroup)

                        # merge the entity ids
                        atom.SetIntData('entity_serial', entity_id_ref)

                        # merge pdb chain asu information
                        atom.SetStringData('pdb_chain_asu_id', pdb_chain_asu_id_ref)

                    app.log.debug("alternate location {0} set for old location {1} "
                                  "of residue {2}".format(altloc_new, altloc, residue))

            app.log.debug("PDB Chain ID set to {0} for alternate location group"
                          .format(pdb_chain_id_ref))

        # add cleaned ligand entity ids to original structure again
        structure.SetData('ligand_entity_serials', list(entity_serials))
