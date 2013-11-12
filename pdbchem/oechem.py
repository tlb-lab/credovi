"""
Used to monkey-patch OpenEye classes.
"""
import os
from functools import partial, total_ordering, wraps
from operator import itemgetter
from collections import Counter # Python 2.7+

from openeye.oechem import *
from openeye.oeiupac import *
from eyesopen.predicates import *

# OEChem input flavour to use for reading PDB entries
# this flavour will keep alternate location groups
OE_PDB_INPUT_FLAVOR = (OEIFlavor_Generic_Default | OEIFlavor_PDB_ALTLOC |
                       OEIFlavor_PDB_Default | OEIFlavor_PDB_DATA)

def _reconstruct_mol(cls, state):
    """
    This reconstructor is used to create a new instance of an OEChem OEMolBase
    that was previously serialised as OEB string.
    """
    old = mol_from_string(state, oeformat=OEFormat_OEB)
    return cls.__call__(old)

def GetListData(self, tag):
    """
    """
    return [data.GetData() for data in self.GetDataIter(OEGetTag(tag))]

def SetListData(self, tag, iterable):
    """
    """
    self.DeleteData(tag)
    for item in iterable: self.AddData(tag, item)

OEGraphMol.__repr__          = lambda self: '<{}({})>'.format(self.__class__.__name__, id(self))
OEGraphMol.__eq__            = lambda self, other: hash(self) == hash(other)
OEGraphMol.__ne__            = lambda self, other: not self == other
OEGraphMol.__lt__            = lambda self, other: OECount(self,OEIsHeavy()) < OECount(other,OEIsHeavy())
OEGraphMol.__le__            = lambda self, other: OECount(self,OEIsHeavy()) <= OECount(other,OEIsHeavy())
OEGraphMol.__gt__            = lambda self, other: OECount(self,OEIsHeavy()) > OECount(other,OEIsHeavy())
OEGraphMol.__ge__            = lambda self, other: OECount(self,OEIsHeavy()) >= OECount(other,OEIsHeavy())
OEGraphMol.__hash__          = lambda self: hash(OECreateIsoSmiString(self))

# Required for pickling. The pickling strategy for OEGraphMolBase objects is simply
# to serialise them as OpenEye binary strings (OEB).
OEGraphMol.__getstate__      = lambda self: mol_to_oeb(self)

# Required for pickling. Reconstructs the OEGraphMolBase object by reading in an OEB
# string (__getstate__) and using it as input for the reconstruct function
# (together with the class to make sure it gets deserialised to it's proper
# original class).
OEGraphMol.__reduce__        = lambda self: (_reconstruct_mol, (self.__class__,
                                                               self.__getstate__()))

OEGraphMol.__getitem__       = lambda self, idx: self.GetAtom(OEHasAtomIdx(idx))
OEGraphMol.__iter__          = lambda self: self.GetAtoms()
OEGraphMol.GetListData       = GetListData
OEGraphMol.SetListData       = SetListData

OEMol.__repr__          = lambda self: '<{}({})>'.format(self.__class__.__name__, id(self))
OEMol.__eq__            = lambda self,other: hash(self) == hash(other)
OEMol.__ne__            = lambda self,other: not self == other
OEMol.__lt__            = lambda self,other: OECount(self,OEIsHeavy()) < OECount(other,OEIsHeavy())
OEMol.__le__            = lambda self,other: OECount(self,OEIsHeavy()) <= OECount(other,OEIsHeavy())
OEMol.__gt__            = lambda self,other: OECount(self,OEIsHeavy()) > OECount(other,OEIsHeavy())
OEMol.__ge__            = lambda self,other: OECount(self,OEIsHeavy()) >= OECount(other,OEIsHeavy())
OEMol.__hash__          = lambda self: hash(OECreateIsoSmiString(self))

# Required for pickling. The pickling strategy for OEMolBase objects is simply
# to serialise them as OpenEye binary strings (OEB).
OEMol.__getstate__      = lambda self: mol_to_oeb(self)

# Required for pickling. Reconstructs the OEMolBase object by reading in an OEB
# string (__getstate__) and using it as input for the reconstruct function
# (together with the class to make sure it gets deserialised to it's proper
# original class).
OEMol.__reduce__        = lambda self: (_reconstruct_mol, (self.__class__,
                                                           self.__getstate__()))

OEMol.__getitem__       = lambda self,idx: self.GetAtom(OEHasAtomIdx(idx))
OEMol.__iter__          = lambda self: self.GetAtoms()
OEMol.GetListData       = GetListData
OEMol.SetListData       = SetListData

OEAtomBase.__repr__     = lambda self: '<OEAtomBase({0}:{1})>'.format(self.GetIdx(), self.GetName())
OEAtomBase.__eq__       = lambda self,other: self.GetIdx() == other.GetIdx()
OEAtomBase.__ne__       = lambda self,other: self.GetIdx() != other.GetIdx()
OEAtomBase.__lt__       = lambda self,other: self.GetIdx() < other.GetIdx()
OEAtomBase.__le__       = lambda self,other: self.GetIdx() <= other.GetIdx()
OEAtomBase.__gt__       = lambda self,other: self.GetIdx() > other.GetIdx()
OEAtomBase.__ge__       = lambda self,other: self.GetIdx() >= other.GetIdx()
OEAtomBase.__hash__     = lambda self: hash(self.GetIdx())

OEResidue.__repr__      = lambda self: '<OEResidue({0[0]} {0[1]} {0[2]}{0[3]})>'.format(self.GetPDBId())
OEResidue.GetPDBId      = lambda self: (self.GetChainID(), self.GetName().strip(), self.GetResidueNumber(), self.GetInsertCode())
OEResidue.__eq__        = lambda self,other: OESameResidue(self,other)
OEResidue.__ne__        = lambda self,other: not OESameResidue(self,other)
OEResidue.__lt__        = lambda self,other: self.GetPDBId() < other.GetPDBId()
OEResidue.__le__        = lambda self,other: self.GetPDBId() <= other.GetPDBId()
OEResidue.__gt__        = lambda self,other: self.GetPDBId() > other.GetPDBId()
OEResidue.__ge__        = lambda self,other: self.GetPDBId() <= other.GetPDBId()
OEResidue.__hash__      = lambda self: hash(self.GetPDBId())

def OEClearNonChiralAtomStereo(molecule):
    """
    Removes false atom stereo information from non-chiral atoms.
    """
    predicate = OEAndAtom(OEHasAtomStereoSpecified(), OENotAtom(OEIsChiralAtom()))

    for atom in molecule.GetAtoms(predicate):
        atom.SetStereo(list(atom.GetAtoms()), OEAtomStereo_Tetrahedral,
                       OEAtomStereo_Undefined)

def OEClearNonChiralBondStereo(molecule):
    """
    Removes false bond stereo information from non-chiral bonds.
    """
    predicate = OEAndBond(OEHasOrder(2), OEAndBond(OEHasBondStereoSpecified(),
                                                   OENotBond(OEIsChiralBond())))

    for bond in molecule.GetBonds(predicate):
        bond.SetStereo([bond.GetBgn(), bond.GetEnd()], OEBondStereo_CisTrans,
                       OEBondStereo_Undefined)

def mol_from_file(path, fmt=None, **kwargs):
    """
    """
    if not os.path.exists(path):
        raise IOError("cannot open file: {} does not exist.".format(path))

    if not os.path.isfile(path):
        raise IOError("cannot open file: {} is not valid.".format(path))

    if kwargs.get('oemol'): method = 'GetOEMols'
    else: method = 'GetOEGraphMols'

    if not fmt: _, fmt = path.rsplit('.', 1)

    ifs = oemolistream()
    if kwargs.get('gz'): ifs.Setgz(gz) # true if input file is gzipped

    if fmt == 'pdb':
        ifs.SetFlavor(OEFormat_PDB, OE_PDB_INPUT_FLAVOR)

    ifs.open(str(path)) # no unicode

    # file contains only a single molecule
    if kwargs.get('multi'): return getattr(ifs, method)()

    # by default, NMR models are treated as sequential molecules, which means
    # only the first model is considered here
    try: molecule = getattr(ifs, method)().next()
    except StopIteration: return None
    else: return molecule

def mol_from_string(string, oeformat=None, **kwargs):
    """
    Returns a new OEMol created from the molecular structure given as a string
    in the specified format. An OEMol object is returned because the molecule
    string might contain a molecule with conformers.

    :param oeformat: File format flag from the OEFormat namespace.
    """
    # raise error if no file format was specified
    if not oeformat: ValueError("cannot read molecule: file format missing.")

    iss = oemolistream()
    iss.SetFormat(oeformat)

    # use the eyesopen default PDB input flavour to keep disordered atoms
    if oeformat == OEFormat_PDB:
        iss.SetFlavor(OEFormat_PDB, OE_PDB_INPUT_FLAVOR)

    if iss.openstring(string):
        try:
            molecule = iss.GetOEMols().next()
            molecule.SetTitle(str(kwargs.get('title','')))
        except StopIteration:
            raise ValueError("string does not contain a valid molecule: {}"
                             .format(string))
        else:
            return molecule

def mol_from_smiles(smiles, canon=False, strict=True):
    """
    Returns an OEGraphMol for the given SMILES string.

    :param canon: Can be used to circumvent the post-processing kekulization test.
                  Passing a Boolean true value to this argument indicates to the
                  parser that the SMILES string should be assumed to be well-formed
                  and the usual kekulization (by calling OEKekulize) step may be
                  omitted. This can be used to speed-up parsing of a large database,
                  but has the side-effect that bond orders are not correctly assigned
                  for aromatic molecules.
    :param strict: Controls whether the parser should operate in strict mode. By
                   default, the SMILES parser attempts to process any reasonably
                   formed SMILES string. If true, the parser applies more rigorous
                   sanity checking.
    """
    if not isinstance(smiles, str):
        raise TypeError("cannot parse SMILES: input must be an instance of string "
                        "but is of type {}".format(type(smiles)))

    molecule = OEGraphMol()

    if OEParseSmiles(molecule, smiles, canon, strict):

        # input SMILES might be kekulized - important to perceive aromaticity here
        OEAssignAromaticFlags(molecule)

        return molecule

def mol_from_iupac_name(name, language=OELanguage_AMERICAN, strict=True):
    """
    """
    molecule = OEGraphMol()
    
    name = str(name.strip())
    name = OEReorderIndexName(name) or name
    
    if not OEParseIUPACName(molecule, name) and strict:
        raise TypeError("cannot parse IUPAC name: {}.".format(name))
        
    return molecule

def mol_to_iupac_name(molecule, style=OENamStyleOpenEye, strict=True):
    """
    """
    name = OECreateIUPACName(molecule, style)

    if strict and 'BLAH' in name:
        raise TypeError("cannot create IUPAC name.")

    return name

def standardise_smiles(smiles, isomeric=True, reset_charges=True, kekule=False):
    """
    Returns the canonical SMILES string of a molecule.
    """
    smiflag = OESMILESFlag_Canonical

    molecule = mol_from_smiles(smiles)

    # reset formal charges
    if reset_charges:
        for atom in molecule.GetAtoms(OEIsPolar()):
            atom.SetImplicitHCount(0)
            atom.SetFormalCharge(0)

        OEAssignImplicitHydrogens(molecule)
        OEAssignFormalCharges(molecule)

    OEFindRingAtomsAndBonds(molecule)
    OEAssignAromaticFlags(molecule,OEAroModelOpenEye)

    # create SMILES with stereo information
    if isomeric:
        smiflag |= OESMILESFlag_ISOMERIC

        OEPerceiveChiral(molecule)
        OEClearNonChiralAtomStereo(molecule)
        OEClearNonChiralBondStereo(molecule)

    if kekule:
        for bond in molecule.GetBonds(OEIsAromaticBond()):
            bond.SetIntType(5)

        OECanonicalOrderAtoms(molecule)
        OECanonicalOrderBonds(molecule)
        OEClearAromaticFlags(molecule)
        OEKekulize(molecule)

    return OECreateSmiString(molecule,smiflag)

def mol_to_smiles(molecule, isomeric=True, reset_charges=True, kekule=False,
                  from3d=False):
    """
    Returns the canonical SMILES string of a molecule.
    """
    smiflag = OESMILESFlag_Canonical

    # reset formal charges
    if reset_charges:

        # hydrogens need to be suppressed
        OESuppressHydrogens(molecule)

        for atom in molecule.GetAtoms(OEIsPolar()):
                atom.SetImplicitHCount(0)
                atom.SetFormalCharge(0)

        OEAssignImplicitHydrogens(molecule)
        OEAssignFormalCharges(molecule)

    OEFindRingAtomsAndBonds(molecule)
    OEAssignAromaticFlags(molecule,OEAroModelOpenEye)

    # get stereo information from 3D coordinates
    if from3d: OE3DToInternalStereo(molecule)

    # create smiles with stereo information
    if isomeric:
        smiflag |= OESMILESFlag_ISOMERIC

        OEPerceiveChiral(molecule)
        OEClearNonChiralAtomStereo(molecule)
        OEClearNonChiralBondStereo(molecule)

    if kekule:
        for bond in molecule.GetBonds(OEIsAromaticBond()):
            bond.SetIntType(5)

        OECanonicalOrderAtoms(molecule)
        OECanonicalOrderBonds(molecule)
        OEClearAromaticFlags(molecule)
        OEKekulize(molecule)

    return OECreateSmiString(molecule,smiflag)

def mol_to_string(molecule, oeformat=None):
    """
    Returns a string containing the structure of the input molecule in the given
    format from the OEFormat namespace.
    """
    oss = oemolostream()
    oss.SetFormat(oeformat)

    oss.openstring()
    OEWriteMolecule(oss, molecule)

    return oss.GetString()

mol_to_pdb = partial(mol_to_string, oeformat=OEFormat_PDB)
mol_to_sdf = partial(mol_to_string, oeformat=OEFormat_SDF)
mol_to_oeb = partial(mol_to_string, oeformat=OEFormat_OEB)

def kekulize_mol(molecule):
    """
    """
    for bond in molecule.GetBonds(OEIsAromaticBond()):
        bond.SetIntType(5)

    OECanonicalOrderAtoms(molecule)
    OECanonicalOrderBonds(molecule)
    OEClearAromaticFlags(molecule)
    OEKekulize(molecule)

    return molecule

def get_subsetmol(molecule, predicate, remove=False):
    """
    Returns the subset of a molecule defined by a given predicate as new molecule.
    """
    subsetmol = OEMol()

    subsetatoms = molecule.GetAtoms(predicate)
    member = OEIsAtomMember(subsetatoms)

    OESubsetMol(subsetmol, molecule, member, False, False)

    # remove the selection from the original molecule
    if remove:
        for atom in molecule.GetAtoms(predicate):
            molecule.DeleteAtom(atom)

    return subsetmol

def reset_charges(molecule):
    """
    Resets the formal charges of a molecule.
    """
    for atom in molecule.GetAtoms(OEIsPolar()):
        atom.SetImplicitHCount(0)
        atom.SetFormalCharge(0)

    OEAssignImplicitHydrogens(molecule)
    OEAssignFormalCharges(molecule)
    OEAssignAromaticFlags(molecule)

    return molecule

def has_undef_stereo(molecule):
    """
    Returns True if the molecule has undefined stereo centers.
    """
    OEPerceiveChiral(molecule)

    predicate = OEAndAtom(OENotAtom(OEHasAtomStereoSpecified()), OEIsChiralAtom())

    return OECount(molecule, predicate) > 0

def has_std_atoms(molecule):
    """
    """
    return OECount(molecule, OEIsStdAtom()) == molecule.NumAtoms()

def disconnected_components(molecule, min_atoms=5, max_atoms=70, remove=False,
                            adjust_hcount=False, rgroup=False, atom_map=False):
    """
    Returns a list of molecules that are disconnected in the input molecule.
    """
    atomarray = OEAtomArray(molecule.GetMaxAtomIdx())
    numparts, partlist = OEDetermineComponents(molecule)
    pred = OEPartPredAtom(partlist)

    components = []
    for i in range(1, numparts+1):
        pred.SelectPart(i)
        component = OEGraphMol()
        OESubsetMol(component, molecule, pred, adjust_hcount, rgroup, atomarray)

        if min_atoms <= component.NumAtoms() <= max_atoms:
            if atom_map: components.append((component, OEAtomArray(atomarray)))
            else: components.append(component)

            # remove the component from the original molecule
            if remove:
                for atom in molecule.GetAtoms(pred):
                    molecule.DeleteAtom(atom)

    return components

def murcko_scaffold(molecule, generic=False):
    """
    """
    OESuppressHydrogens(molecule);

    librarygen = OELibraryGen();
    librarygen.Init("[*;R:1]-[A;!#1:2]>>[*:1].[*:2]");
    librarygen.SetExplicitHydrogens(False);
    librarygen.SetValenceCorrection(False);

    hits = librarygen.SetStartingMaterial(molecule, 0);

    # if there are no hits then molecule is already murcko scaffold
    if hits > 0:

        # iterate through all reaction products
        for product in librarygen.GetProducts():

            # cleave the reaction product to determine if fragment is part of scaffold or not
            for fragment, atom_map in disconnected_components(product, min_atoms=0, atom_map=True):

                # fragment does not contain any rings and can be removed
                if OECount(fragment, OEAtomIsInRing()) == 0:

                    # remove functional groups from original molecule
                    for idx, fragatom in enumerate(atom_map):
                        if fragatom:
                            atom = molecule.GetAtom(OEHasAtomIdx(idx))
                            molecule.DeleteAtom(atom)

        # make a generic scaffold: only carbon atoms and single bonds
        if generic:
            for atom in molecule.GetAtoms():
                atom.SetAtomicNum(OEElemNo_C);
                atom.SetAromatic(False);

            # set single bonds and remove aromaticity information
            for bond in molecule.GetBonds():
                bond.SetOrder(1);
                bond.SetAromatic(False);

    # normalise molecule
    reset_charges(molecule)

    # create smiles representation of the new murcko scaffold
    return OECreateIsoSmiString(molecule);

def mcs(m1, m2, min_atoms=5, max_matches=64, unique=True,
        atomexpr=OEExprOpts_DefaultAtoms, bondexpr=OEExprOpts_DefaultBonds,
        mcstype=OEMCSType_Approximate, scorefunc='max_atoms_complete_cycles'):
    """
    Returns a list of Maximum Common Substructures (MCS) shared between the two
    molecules.
    """
    # create maximum common substructure object with default atoms and bondS
    mcss = OEMCSSearch(m1, atomexpr, bondexpr, mcstype)

    # scoring function to use for MCS search
    if scorefunc=='max_atoms':
        scorefunc = OEMCSMaxAtoms()
    elif scorefunc=='max_bonds':
        scorefunc = OEMCSMaxBonds()
    elif scorefunc=='max_atoms_complete_cycles':
        scorefunc = OEMCSMaxAtomsCompleteCycles()
    elif scorefunc=='max_bonds_complete_cycles':
        scorefunc = OEMCSMaxBondsCompleteCycles()
    else:
        raise ValueError("{} is not an allowed option for scorefunc."
                         .format(scorefunc))

    # set scoring function
    mcss.SetMCSFunc(scorefunc)

    # ignore matches smaller than the minimum number of heavy atoms
    mcss.SetMinAtoms(min_atoms)
    mcss.SetMaxMatches(max_matches)

    hits = set()
    for match in mcss.Match(m2, unique):
        substructure = OEGraphMol()
        OESubsetMol(substructure, match, True)
        hits.add(substructure)

    return hits

def multimcs(data, threshold=1.0, min_atoms=5, atomexpr=16843264,
             bondexpr=OEExprOpts_DefaultBonds, mcstype=OEMCSType_Approximate,
             scorefunc='max_atoms_complete_cycles', kekulize=True):
    """
    Determines the Maximum Common Substructure (MCS) in a set of molecules in
    SMILES format.

    :param smistrings: list (or iterable) of molecules in SMILES format.
    :param threshold: minimum number of compounds that must share the MCS: defaults
                      to 1.0, which means all. 0.8 For example would mean that
                      at least 80% of compounds must share the MCS.
    :param min_atoms: minimum number of atoms the MCS must have.
    :param atomexpr: the OpenEye MCS atom expression to use: defaults to
                     OEExprOpts_AtomicNumber | OEExprOpts_Aromaticity |
                     OEExprOpts_EqCHalogen.
    :param bondexpr: the OpenEye MCS bond expression to use: defaults to
                     OEExprOpts_DefaultBonds.
    :param mcstype: OpenEye MCS type to use, defaults to OEMCSType_Approximate,
                    i.e. an approximate search is used. Should be sufficient for
                    almost all cases.
    :param scorefunc: OpenEye MCS scoring function to use; defaults to
                      OEMCSMaxAtomsCompleteCycles.

    ..important:: Requires Python 2.7+.
    """
    # create molecule objects (might consume a lot of memory)
    compounds = []

    # get the molecules from the file if provided
    if isinstance(data, str) and os.path.isfile(data):
        ifs = oemolistream()
        ifs.open(data)

        # add molecules to temporary list of molecules
        for molecule in ifs.GetOEGraphMols():
            compounds.append(OEGraphMol(molecule))

    # the input data is a list of SMILES strings
    elif isinstance(data, list):
        compound = OEGraphMol()

        # convert each SMILES string to a molecule and add to list
        for smiles in data:
            if OEParseSmiles(compound, smiles):
                compounds.append(OEGraphMol(compound))
                compound.Clear()

            else:
                raise ValueError("cannot parse SMILES: {}".format(smiles))

    # return the original molecule if only one is given
    if len(compounds) == 1:
        return OECreateSmiString(compounds[0])

    # create a new max() function that can return 0 on an empty list
    maximum = partial(reduce, max)

    # freeze the MCS method with the specified parameters
    mcsfunc = partial(mcs, min_atoms=min_atoms, atomexpr=atomexpr,
                      bondexpr=bondexpr, mcstype=mcstype, scorefunc=scorefunc)

    # substructure search object
    subsearch = OESubSearch()

    # sort compounds by size
    # the largest will be the first query for MCS searching
    compounds.sort(key=lambda molecule: molecule.NumAtoms())

    # determine the minimum number of expected matches: defaults to all compounds
    min_exp_matches = len(compounds) * threshold

    # initialize a pool of possible MCS's with the largest compound as reference
    substructures = set([compounds.pop()])

    # create a counter for MCS's: the one we want will be equal to the number of
    # compounds, i.e. found in all structures
    matchcount = Counter()

    # shortcut to the values of the match counter
    hits = matchcount.viewvalues()

    # loop until an MCS is found in the minimum number of expected molecules
    # we need to do +1 here because we removed a compound to serve as initial query
    while maximum(hits, 0) + 1 < min_exp_matches:

        # end immediately if no substructure (and possible MCS) is left in the pool
        if not substructures: return None

        # get a substructure match from the pool as a reference
        substruct = substructures.pop()

        # initialize substructure search for picked reference
        subsearch.Init(substruct, atomexpr, bondexpr)

        # iterate through compounds and record number of hits
        for compound in compounds:

            # check for a simple match first
            if subsearch.SingleMatch(compound):
                matchcount[substruct] += 1

                # break here if an MCS is now found in the minimum number of
                # compounds in the input set
                if maximum(hits, 0) + 1 >= min_exp_matches: break

            # otherwise get the MCS between the reference and the current target
            # and add the new hits to our list of substructures as references
            # for the next search
            else:
                for hit in mcsfunc(substruct, compound):
                    substructures.add(hit)

                # cancel search if one already fails to safe time: the MCS by
                # default has to match every compound, so there is no point to
                # test the others if this does not contain this substructure
                # unless a threshold less than 100% was given
                if not threshold < 1.0: break

    # take the first substructure of the hit with the largest number of atoms
    result = matchcount.most_common(1)

    if result:
        molecule, count = result.pop()

        if kekulize: result = kekulize_mol(molecule)
        return OECreateIsoSmiString(molecule)
