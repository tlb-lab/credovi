import os
import logging
from itertools import groupby, imap

from openeye.oespicoli import OESurface, OEMakeMolecularSurface
from eyesopen.oechem import *

from credovi import app # to get configuration and logger
from credovi.util.timer import timer

def get_structure(path, pdbcode=None):
    """
    """
    if os.path.exists(path):
        ifs = oemolistream(path)

        if ifs.IsValid():
            structure = OEGraphMol()
            success = OEReadMolecule(ifs, structure)

            if success:
                if pdbcode: structure.SetTitle(str(pdbcode))
                return structure

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
            assembly_dir = os.path.join(app.config.get('directories','quat_oeb'),
                                        pdb[1:3].lower(), pdb.lower())

            if os.path.exists(assembly_dir):

                # add all assemblies in directory to our assembly set
                for oeb in os.listdir(assembly_dir):
                    assemblyset.append(oeb)

                # add assembly set to all the other assembly sets
                if assemblyset:
                    assemblysets.append((pdb.upper(), assemblyset))

            # no assemblies exists for the given PDB code
            else:
                app.log.warn("no assembly set found for PDB entry: {pdb}"
                             .format(pdb=pdb))

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
            pdbs = map(str.upper, app.config['test sets'][args.testset])

        assemblysets = _get_assembly_sets(pdbs)

    # get only the asymmetric unit structures
    else:
        assemblysets = []

    return sorted(assemblysets)

def get_aromatic_rings(structure, patterns):
    """
    """
    hit = 0
    subsearch = OESubSearch()

    for pattern in patterns:
        subsearch.Init(str(pattern))
        subsearch.SetMaxMatches(0)

        for match in subsearch.Match(structure, True): # only unique matches
            hit += 1

            # get all the atoms that are part of this aromatic ring
            ringatoms = list(match.GetTargetAtoms())

            # get the residue information for each atom to make sure that they
            # all belong to the same. Aromatic rings with more than one residue
            # cannot be stored in CREDO.
            residues = {OEAtomGetResidue(atom) for atom in ringatoms}

            # log an appropriate message and skip this ring
            if len(residues) > 1:
                app.log.warn("cannot add aromatic ring {}: belongs to more than "
                             "one residue ({}).".format(hit, ', '.join(str(res) for res in residues)))
                continue

            yield (hit, ringatoms)


def get_pi_groups(structure, patterns):
    """
    """
    hit = 0
    subsearch = OESubSearch()

    for pattern in patterns:
        subsearch.Init(str(pattern))
        subsearch.SetMaxMatches(0)

        for match in subsearch.Match(structure, True): # only unique matches
            hit += 1

            # get all the atoms that are part of this pi system
            pi_atoms = list(match.GetTargetAtoms())

            yield (hit, pi_atoms)


def get_ligands(structure):
    """
    Returns a list of ligand molecules found in the structure.
    """
    ligands = []

    # check if structure has any ligands at all
    app.log.debug("checking for ligands in structure {}...".format(structure.GetTitle()))
    if structure.HasData('ligand_entity_serials'):
        app.log.debug("structure {} has ligands. processing them.".format(structure.GetTitle()))

        ### create a mapping between atom ids and ligand entity serials for partitioning

        entity_serials = structure.GetData('ligand_entity_serials')
        partlist = OEIntArray(structure.GetMaxAtomIdx())

        # create ligand partition list
        for atom in structure.GetAtoms():
            entity_serial = atom.GetIntData('entity_serial')

            if entity_serial in entity_serials: partlist[atom.GetIdx()] = entity_serial

        ligpred = OEPartPredAtom(partlist)

        # create a new molecule for each ligand in the structure
        for entity_serial in entity_serials:
            ligpred.SelectPart(entity_serial)

            ligand = OEGraphMol()
            OESubsetMol(ligand, structure, ligpred)

            # remove pdb header from ligand
            OEClearPDBData(ligand)
            OEDeleteEverythingExceptTheFirstLargestComponent(ligand)

            # add the ligand to list
            if ligand.NumAtoms():

                # get all the residues (components) of this ligand
                residues = set((OEAtomGetResidue(atom) for atom in ligand.GetAtoms()))

                # sort residues by residue number
                residues = sorted(residues, key=lambda r: r.GetResidueNumber())

                # add information for single component ligand
                if len(residues) == 1:
                    residue = residues.pop()
                    pdb_chain_id, name, res_num = residue.GetPDBId()[:-1]

                    ligand.SetIntData('res_num', res_num)

                # add information for peptide ligand only
                elif len(residues) <= 10:
                    pdb_chain_id = residues[0].GetChainID()
                    name = '-'.join((residue.GetName().strip() for residue in residues))

                    app.log.debug("chain {0} is peptide ligand with sequence {1}."
                                  .format(pdb_chain_id, name))

                # ligand has too many residues - must be an error
                else:
                    app.log.fatal("PDB {0}: ligand with entity serial number {1} has {2} residues - aborting!:\n{3}"
                                  .format(structure.GetTitle(), entity_serial, len(residues), residues))
                    app.close()

                # add required information to ligand molecule
                ligand.SetTitle(name)
                ligand.SetStringData('name', name)
                ligand.SetStringData('pdb_chain_id', pdb_chain_id)
                ligand.SetIntData('entity_serial', entity_serial)

                ligands.append(ligand)

            # something must have gotten wrong with the entity serial number
            else:
                app.log.warn('ligand with entity serial number {0} cannot be found in structure!'
                             .format(entity_serial))

    return ligands

@timer("all atom types set in {0:.2f} seconds.")
def set_atom_type_flags(structure, set_ambig=False):
    """
    Sets the CREDO atom type flags to all atoms in the given structure. Make sure
    that all atom types are set in lowercase.
    """
    atom_types = app.config['atom types']

    pdb_atom_types = app.config['pdb_atom_types']
    pdb_atom_types_ext = app.config['pdb_atom_types_ext']

    for atom in structure.GetAtoms(OEIsStdProteinResidue(strict=True)):
        atomres = OEAtomGetResidue(atom)

        atom_key = "{resname}_{atomname}".format(
            resname=atomres.GetName().strip(),
            atomname=atom.GetName().strip())

        pdb_type = pdb_atom_types.get(atom_key, pdb_atom_types_ext.get(atom_key) if set_ambig else None)

        if pdb_type:
            if isinstance(pdb_type, basestring):
                atom.SetIntData(str(pdb_type), 1)
            else:
                for typ in pdb_type:
                    atom.SetIntData(str(typ), 1)
            atom.SetIntData('prot_atom', 1)


    for atom_type, smartsdict in atom_types.items():

        # only for SMARTS patterns that hit single atoms
        if atom_type not in ('pos ionisable','neg ionisable','pi-system') and smartsdict:
            for smarts in smartsdict.values():
                for atom in structure.GetAtoms(OEAndAtom(OENotAtom(OEIsStdProteinResidue(strict=True)),
                                                         OEMatchAtom(str(smarts)))):
                    atom.SetIntData(str(atom_type),1)
                    if atom_type == "hbond acceptor":
                        atom.SetIntData("xbond acceptor",1)

    # some of the following smarts hit multiatom groups therefore pattern matching has to be used
    ss = OESubSearch()

    # positively ionisable atoms
    for smarts in atom_types['pos ionisable'].values():
        ss.Init(str(smarts))
        ss.SetMaxMatches(0)

        for match in ss.Match(structure):
            for atom in match.GetTargetAtoms():
                if atom.HasData('prot_atom'):
                    continue
                atom.SetIntData('pos ionisable',1)

    for atom in structure.GetAtoms(OEIsMetal()):
        atom.SetIntData('pos ionisable',1)

    # negatively ionisable atoms
    for smarts in atom_types['neg ionisable'].values():
        ss.Init(str(smarts))
        ss.SetMaxMatches(0)

        for match in ss.Match(structure):
            for atom in match.GetTargetAtoms():
                if atom.HasData('prot_atom'):
                    continue
                atom.SetIntData('neg ionisable',1)

    # all water molecules are hydrogen bond donors and acceptors
    for atom in structure.GetAtoms(OEIsWater()):
        atom.SetIntData('hbond acceptor',1)
        atom.SetIntData('hbond donor',1)

@timer("disconnected components identified in {0:.2f} seconds.")
def identify_disconnected_components(structure):
    """
    Identifies all disconnected components in this structure and assigns an entity
    serial number to each. The serial number will be used later on to identify
    intramolecular interactions.
    """
    numparts, partlist = OEDetermineComponents(structure)
    pred = OEPartPredAtom(partlist)

    for atom in structure.GetAtoms(OENotIsWater):
        atom.SetIntData('entity_serial', partlist[atom.GetIdx()])

    app.log.debug("structure contains {0} disconnected components.".format(numparts))

    return numparts, pred

@timer("CA/CB ratio calculated in {0:.2f} seconds.")
def calc_calpha_ratio(structure):
    """
    Returns True if any of the polypeptide chains has a CB/CA ratio below 0.8.
    """
    CALPHA_FLAG = False

    # create predicate to only include polypeptides (without glycine)
    functors = (OEIsStdProteinResidue(), OENotIsWater, OENotAtom(OEIsInGLY), OEIsHeavy())
    predicate = reduce(OEAndAtom, functors)

    atoms = structure.GetAtoms(predicate)

    # sort atoms by chain, otherwise groupby will produce false results
    atoms = sorted(atoms, key=lambda atom: OEAtomGetResidue(atom).GetChainID())

    for pdb_chain_id, atoms in groupby(atoms, lambda atom: OEAtomGetResidue(atom).GetChainID()):
        atoms = list(atoms)

        # get calpha and cbeta count
        Ca = len(filter(OEHasPDBAtomIndex(OEPDBAtomName_CA), atoms))
        Cb = len(filter(OEHasPDBAtomIndex(OEPDBAtomName_CB), atoms))

        # check if chain contains ca atoms
        try:
            ratio = float(Cb) / float(Ca)
        except ZeroDivisionError:
            app.log.warn('chain {0} does not contain any CA atoms.'.format(pdb_chain_id))
            continue

        # debug cb/ca ratio info
        app.log.debug("CB/CA ratio for chain {0}: {1:.2f}.".format(pdb_chain_id,ratio))

        # get the number of residues in the chain
        residues = list(set(imap(lambda x: OEAtomGetResidue(x),atoms)))

        # treat short chain as ligand
        if len(residues) > 11 and Ca and ratio < 0.8: CALPHA_FLAG = True

    return CALPHA_FLAG

@timer("entity serial numbers assigned to all atoms in ASU in {0:.2f} seconds.")
def assign_entity_serials(structure, pdb_polymer_info, pdb_ligand_info, pdb_sstruct_info):
    """
    Assigns a serial number to each identified entity in the structure. An entity
    is one (or more) residue(s) that form a distinct biological molecule, e.g. a
    polypeptide chain, ligand, nucleic acid, etc.

    This function also assigns a bitmask (entity_type_bm) in the form below to each
    atom:

        Entity bitmask
            is_solvent      1
            is_ligand       2
            is_saccharide   4
            is_rna          8
            is_dna          16
            is_protein      32
    """
    ligand_entity_serials = []

    ### create numeric mapping for all entities in structure

    # create entity serial for all pdb chains
    pdb_chain_ids = sorted(set(row[0] for row in pdb_polymer_info.keys()))
    entities = dict(((pdb_chain_id, None, None, ' '), entity_serial) for entity_serial, pdb_chain_id in enumerate(pdb_chain_ids,1))

    for entry, entity_serial in entities.items():
        app.log.debug('entity serial number {0} assigned to polymer chain {1}.'
                      .format(entity_serial, entry[0]))

    # get the current highest entity serial number
    if entities: max_entity_serial = max(entities.values())

    # pdb entry does not contain any polymer chains (1ao4)
    else: max_entity_serial = 0

    # create entity serial for all ligands
    for entry in pdb_ligand_info:

        # ligand is whole chain / peptide ligand
        if entry in entities:

            # store peptide ligand entity serial number info
            ligand_entity_serials.append(entities[entry])

            # debug message
            app.log.debug('entity serial number {0} now classified as peptide ligand.'
                         .format(entities[entry]))

        # ligand is single residue
        else:

            # increment entity serial number for ligand
            max_entity_serial += 1

            # update entity mapping
            entities[entry] = max_entity_serial

            # store ligand entity serial number info
            ligand_entity_serials.append(max_entity_serial)

            # debug message
            app.log.debug('entity serial number {0} assigned to ligand {1}'
                          .format(max_entity_serial, ' '.join(map(str,entry))))

    ### assign entity serial number, bitmask and secondary structure info to each atom

    # create a combined residue mapping
    pdb_res_info = dict(pdb_polymer_info, **pdb_ligand_info)

    # iterate through all atoms in the structure
    for atom in structure.GetAtoms():

        # get the residue object
        residue = OEAtomGetResidue(atom)
        pdb_id = residue.GetPDBId()

        # add secondary structure serial number
        atom.SetIntData('sstruct_serial', pdb_sstruct_info.get(residue.GetPDBId(),0))

        ## assign entity serial number and entity type bitmask

        # residue belongs to ligand
        if pdb_id in pdb_ligand_info:
            atom.SetIntData('entity_serial', entities[pdb_id])
            atom.SetIntData('entity_type_bm',2)

        # residue is part of a peptide ligand
        elif pdb_id in pdb_polymer_info and (pdb_id[0], None, None, ' ') in pdb_ligand_info:
            atom.SetIntData('entity_serial', entities[(pdb_id[0], None, None, ' ')])
            atom.SetIntData('entity_type_bm',34)

        # residue is part of a polymer
        elif pdb_id in pdb_polymer_info:
            atom.SetIntData('entity_serial', entities[(pdb_id[0], None, None, ' ')])
            atom.SetIntData('entity_type_bm', pdb_res_info.get(pdb_id,0))

        # residue is water (or deuterium)
        elif residue.GetName() in ('HOH','DOD'):
            atom.SetIntData('entity_type_bm',1)

        # should not enter here
        else:
            app.log.warn(("\033[93mno entity serial number and bit mask assigned to atom {0} of {1}. " +
                         "PDB_ID: {2}, size pdb_ligand_info: {3}, size pdb_polymer_info: {4}. Maybe outdated tables?\033[0m")
                         .format(atom, residue, pdb_id, len(pdb_ligand_info), len(pdb_polymer_info)))
            raise RuntimeError("Missing PDB entry on pdb_polymer|ligand_info. Outdated tables?")

    # attach the entity identifiers of ligands (if any) to the structure
    if ligand_entity_serials:
        structure.SetData('ligand_entity_serials', ligand_entity_serials)
        app.log.debug("assigned the following entity serial numbers to all ligands in the ASU: {0}"
                      .format(ligand_entity_serials))

@timer("entity surface atoms were identified in {0:.2f} seconds.")
def identify_surface_atoms(structure):
    """
    Generates the molecular surface of all entities in the structure to identify
    and label all atoms that are exposed to the surface. The resolution is the grid
    spacing to use during surface construction.
    """
    def transform_atom_map(atom_map):
        """
        This generator transforms the boolean atom mapping between a molecule and
        a subset into a tuple list of corresponding atom identifiers.
        """
        entity_atom_idx = -1 # ATOM IDS START AT 0
        for structure_atom_idx, is_entity_atom in enumerate(atom_map):
            if is_entity_atom:
                entity_atom_idx += 1
                yield (entity_atom_idx, structure_atom_idx)

    # boolean mapping that is used to map subset atoms to the original structure
    atom_map = OEAtomArray(structure.GetMaxAtomIdx())

    ### partition the structure into entities

    entity_serials = set()
    partlist = OEIntArray(structure.GetMaxAtomIdx())

    # create the partition map
    for atom in structure.GetAtoms(OENotIsWater):
        entity_serial = atom.GetIntData('entity_serial')
        partlist[atom.GetIdx()] = entity_serial
        entity_serials.add(entity_serial)

    entitypred = OEPartPredAtom(partlist)

    surface = OESurface()

    # create a new molecule for each entity in the structure
    for entity_serial in entity_serials:
        entitypred.SelectPart(entity_serial)

        # create a new molecule for the entity
        entity = OEGraphMol()
        OESubsetMol(entity, structure, entitypred, False, False, atom_map)

        # map entity subsetmol surface atoms back to original structure
        entity_atom_map = dict(transform_atom_map(atom_map))

        # not necessary for single atom entities
        if entity.NumAtoms() > 1:

            # adjust the grid spacing to the size of the entity
            if entity.NumAtoms() > 70: resolution = 1.0
            else: resolution = 0.5

            # make the surface with the given resolution and probe radius
            OEMakeMolecularSurface(surface, entity, resolution, 1.4)

            # map the vertices to their contributing atoms from entity molecule
            vertice_map = OEUIntArray(surface.GetNumVertices())
            surface.GetAtoms(vertice_map)
            surface_atoms = set(vertice_map) # get a non-redundant set of atoms

            # iterate through all surface atoms of the entity, trace them back to the
            # original structure and label them
            for entity_atom_idx in surface_atoms:
                structure_atom_idx = entity_atom_map[entity_atom_idx]
                structure_atom = structure.GetAtom(OEHasAtomIdx(structure_atom_idx))

                # label the entity atom as being exposed to the surface
                structure_atom.SetIntData('is_exposed', 1)

            surface.Clear()

        # label single atom entities as solvent-exposed
        else:
            for atom in entity.GetAtoms():
                structure_atom_idx = entity_atom_map[atom.GetIdx()]
                structure_atom = structure.GetAtom(OEHasAtomIdx(structure_atom_idx))

                # label the entity atom as being exposed to the surface
                structure_atom.SetIntData('is_exposed', 1)

def parse_header(structure):
    """
    Parsed the header of a PDB file.
    """
    MONTH_TO_INT = {"JAN":1,"FEB":2,"MAR":3,"APR":4,"MAY":5, "JUN":6,
                    "JUL":7,"AUG":8,"SEP":9,"OCT":10,"NOV":11,"DEC":12}

    header = {'REVDAT':'','REMARK 350':''}

    # extract information from the PDB header
    for pair in OEGetPDBDataPairs(structure):
        tag, value = pair.GetTag(), pair.GetValue()

        # extract the biomolecule information
        if tag == 'REMARK' and value.startswith(' 350'):
            header['REMARK 350'] += value[5:] + '\n'

        # pdb file is split into several entries
        if tag == 'SPLIT ':
            splits = value.strip().split()
            structure.SetBoolData('is_split', True)

            app.log.warn('split PDB entries: %s' % (splits))

        elif tag == 'REVDAT':
            revision = value[3]
            date = value[7:16]

            # attach first deposition to structure
            if revision == '1':
                day, month, year = date.split('-')
                month = str(MONTH_TO_INT[month])

                if int(year) > 70: year = '19' + year
                else: year = '20' + year

                date = '-'.join([year,month,day])
                structure.SetStringData('deposition', date)

                header['REVDAT'] = date

    return header

def get_ligand_entity_serials(structure):
    """
    Returns a list containing all the entity serial numbers that are assigned to
    the ligands in this structure.
    """
    ligand_entity_serials = set(atom.GetIntData('entity_serial')
                            for atom in structure.GetAtoms(OENotAtom(OEIsWater()))
                            if atom.GetIntData('entity_type_bm') & 2)

    return sorted(ligand_entity_serials)

def extract_entity_with_selection(structure, selection):
    """
    Returns the molecule that matches the selection and extracts it from the
    parent structure. The selection is a string in the form
    <PDB CHAIN ID>.[<RESNUM>] where the residue number is optional.
    """
    pdb_res_num, predicate = '', None
    parts = selection.split('.')

    # try to extract the PDB Chain identifier and residue number
    if len(parts) == 1: pdb_chain_id = parts.pop()
    elif len(parts) == 2: pdb_chain_id, pdb_res_num = parts
    else: raise ValueError("malformed selection: {}.".format(selection))

    # compile the predicate to extract the entity
    if pdb_chain_id.isalpha():
        predicate = OEHasChainID(str(pdb_chain_id))
    if pdb_res_num and pdb_res_num.isdigit():
        predicate = OEAndAtom(predicate, OEHasResidueNumber(int(pdb_res_num)))

    # extract the entity from the parent structure
    entity = get_subsetmol(structure, predicate, remove=True)

    return entity
