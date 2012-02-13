"""
This module, most importantly, monkey-patches the OEGraphMol, OEMol, OEResidue
and OEAtomBase objects. Credovi will not work properly if the patches are missing.
These patches include comparison and hashing methods to make the OpenEye base
classes work with Python sets.
New Predicate functions that are heavily used in Credovi are defined here as well.
"""

from __future__ import absolute_import

import os

from openeye.oechem import *
from openeye.oequacpac import OESetNeutralpHModel
from openeye.oespicoli import (OEMakeAccessibleSurface, OESurface, OESurfaceArea,
                               OESurfaceCliqueArea)

STD_ELEMENTS = frozenset([OEElemNo_C, OEElemNo_N, OEElemNo_O, OEElemNo_S, OEElemNo_P,
                          OEElemNo_Cl, OEElemNo_Br, OEElemNo_I, OEElemNo_F, OEElemNo_Si,
                          OEElemNo_H])

AMINOACIDS = frozenset([OEResidueIndex_ALA, OEResidueIndex_ARG, OEResidueIndex_ASN,
                        OEResidueIndex_ASP, OEResidueIndex_CYS, OEResidueIndex_GLN,
                        OEResidueIndex_GLU, OEResidueIndex_GLY, OEResidueIndex_HIS,
                        OEResidueIndex_ILE, OEResidueIndex_LEU, OEResidueIndex_LYS,
                        OEResidueIndex_MET, OEResidueIndex_PHE, OEResidueIndex_PRO,
                        OEResidueIndex_SER, OEResidueIndex_THR, OEResidueIndex_TRP,
                        OEResidueIndex_TYR, OEResidueIndex_VAL, OEResidueIndex_ASX,
                        OEResidueIndex_GLX, OEResidueIndex_CYX, OEResidueIndex_CYH,
                        OEResidueIndex_HID, OEResidueIndex_HIE, OEResidueIndex_HIP])

NUCLEOTIDES = frozenset([OEResidueIndex_A, OEResidueIndex_C, OEResidueIndex_G,
                         OEResidueIndex_I,OEResidueIndex_T, OEResidueIndex_U,
                         OEResidueIndex_DA, OEResidueIndex_DC, OEResidueIndex_DG,
                         OEResidueIndex_DI,OEResidueIndex_DT, OEResidueIndex_DU])

UNKNOWN_RES = frozenset([OEResidueIndex_UNK, OEResidueIndex_UNL])

# ONLY TEMPORARY!
def GetListData(self, tag):
    """
    """
    return [data.GetData() for data in self.GetDataIter(OEGetTag(tag))]

def SetListData(self, tag, list):
    """
    """
    self.DeleteData(tag)
    for item in list: self.AddData(tag, item)

OEGraphMol.__repr__     = lambda self: '<OEGraphMol({0})>'.format(id(self))
OEGraphMol.__eq__       = lambda self,other: OEExactGraphMatch(self,other)
OEGraphMol.__ne__       = lambda self,other: not OEExactGraphMatch(self,other)
OEGraphMol.__lt__       = lambda self,other: OECount(self,OEIsHeavy()) < OECount(other,OEIsHeavy())
OEGraphMol.__le__       = lambda self,other: OECount(self,OEIsHeavy()) <= OECount(other,OEIsHeavy())
OEGraphMol.__gt__       = lambda self,other: OECount(self,OEIsHeavy()) > OECount(other,OEIsHeavy())
OEGraphMol.__ge__       = lambda self,other: OECount(self,OEIsHeavy()) >= OECount(other,OEIsHeavy())
OEGraphMol.__hash__     = lambda self: hash(mol_to_smiles(self))
OEGraphMol.__getitem__  = lambda self,idx: self.GetAtom(OEHasAtomIdx(idx))
OEGraphMol.__iter__     = lambda self: self.GetAtoms()
OEGraphMol.GetListData  = GetListData # TEMPORARY HACK
OEGraphMol.SetListData  = SetListData # TEMPORARY HACK

OEMol.__repr__          = lambda self: '<OEMol({0})>'.format(id(self))
OEMol.__eq__            = lambda self,other: OEExactGraphMatch(self,other)
OEMol.__ne__            = lambda self,other: not OEExactGraphMatch(self,other)
OEMol.__lt__            = lambda self,other: OECount(self,OEIsHeavy()) < OECount(other,OEIsHeavy())
OEMol.__le__            = lambda self,other: OECount(self,OEIsHeavy()) <= OECount(other,OEIsHeavy())
OEMol.__gt__            = lambda self,other: OECount(self,OEIsHeavy()) > OECount(other,OEIsHeavy())
OEMol.__ge__            = lambda self,other: OECount(self,OEIsHeavy()) >= OECount(other,OEIsHeavy())
OEMol.__hash__          = lambda self: hash(mol_to_smiles(self))
OEMol.__getitem__       = lambda self,idx: self.GetAtom(OEHasAtomIdx(idx))
OEMol.__iter__          = lambda self: self.GetAtoms()
OEMol.GetListData       = GetListData # TEMPORARY HACK
OEMol.SetListData       = SetListData # TEMPORARY HACK

OEResidue.__repr__      = lambda self: '<OEResidue({0[0]} {0[1]} {0[2]}{0[3]})>'.format(self.GetPDBId())
OEResidue.GetPDBId      = lambda self: (self.GetChainID(), self.GetName().strip(), self.GetResidueNumber(), self.GetInsertCode())
OEResidue.__eq__        = lambda self,other: OESameResidue(self,other)
OEResidue.__ne__        = lambda self,other: not OESameResidue(self,other)
OEResidue.__lt__        = lambda self,other: self.GetPDBId() < other.GetPDBId()
OEResidue.__le__        = lambda self,other: self.GetPDBId() <= other.GetPDBId()
OEResidue.__gt__        = lambda self,other: self.GetPDBId() > other.GetPDBId()
OEResidue.__ge__        = lambda self,other: self.GetPDBId() <= other.GetPDBId()
OEResidue.__hash__      = lambda self: hash(self.GetPDBId())

OEAtomBase.__repr__     = lambda self: '<OEAtomBase({0}:{1})>'.format(self.GetIdx(), self.GetAtomicNum())
OEAtomBase.__eq__       = lambda self,other: self.GetIdx() == other.GetIdx()
OEAtomBase.__ne__       = lambda self,other: self.GetIdx() != other.GetIdx()
OEAtomBase.__lt__       = lambda self,other: self.GetIdx() < other.GetIdx()
OEAtomBase.__le__       = lambda self,other: self.GetIdx() <= other.GetIdx()
OEAtomBase.__gt__       = lambda self,other: self.GetIdx() > other.GetIdx()
OEAtomBase.__ge__       = lambda self,other: self.GetIdx() >= other.GetIdx()
OEAtomBase.__hash__     = lambda self: hash(self.GetIdx())


class OEAtomHasIntData(OEUnaryAtomPred):
    '''
    '''
    def __init__(self, pair):
        OEUnaryAtomPred.__init__(self)
        self.pair = pair

    def __call__(self, atom):
        return atom.GetIntData(self.pair[0]) == self.pair[1]

    def CreateCopy(self):
        return OEAtomHasIntData(self.pair).__disown__()

class OEAtomHasAnyIntData(OEUnaryAtomPred):
    '''
    '''
    def __init__(self, data):
        OEUnaryAtomPred.__init__(self)
        self.data = data

    def __call__(self, atom):
        return atom.GetIntData(self.data) > 0

    def CreateCopy(self):
        return OEAtomHasAnyIntData(self.data).__disown__()

class OEAtomBinaryAndIntData(OEUnaryAtomPred):
    """
    Predicate function that can be used to do binary and (&) callbacks on atoms
    having the specified tag.
    """
    def __init__(self, pair):
        OEUnaryAtomPred.__init__(self)
        self.pair = pair

    def __call__(self, atom):
        return atom.GetIntData(self.pair[0]) & self.pair[1]

    def CreateCopy(self):
        return OEAtomBinaryAndIntData(self.pair).__disown__()

class OEAtomBinaryNotIntData(OEUnaryAtomPred):
    '''
    '''
    def __init__(self, pair):
        OEUnaryAtomPred.__init__(self)
        self.pair = pair

    def __call__(self, atom):
        return atom.GetIntData(self.pair[0]) & self.pair[1] == 0

    def CreateCopy(self):
        return OEAtomBinaryNotIntData(self.pair).__disown__()

class OEHasResidueIndex(OEUnaryAtomPred):
    '''
    '''
    def __init__(self, index):
        OEUnaryAtomPred.__init__(self)
        self.index = index

    def __call__(self, atom):
        return OEGetResidueIndex(atom) == self.index

    def CreateCopy(self):
        return OEHasResidueIndex(self.index).__disown__()

class OEIsResidue(OEUnaryAtomPred):
    '''
    Returns True if the atom belongs to the given residue.
    '''
    def __init__(self, residue):
        OEUnaryAtomPred.__init__(self)
        self.residue = residue

    def __call__(self, atom):
        return OESameResidue(self.residue, OEAtomGetResidue(atom))

    def CreateCopy(self):
        return OEIsResidue(self.residue).__disown__()

class OEIsStdAtom(OEUnaryAtomPred):
    '''
    Returns True if the atom is a standard organic element.
    '''
    def __call__(self, atom):
        return atom.GetAtomicNum() in STD_ELEMENTS

class OEIsHetatm(OEUnaryAtomPred):
    '''
    Returns True if the atom is a standard organic element.
    '''
    def __call__(self, atom):
        return OEAtomGetResidue(atom).IsHetAtom()

class OEIsWater(OEUnaryAtomPred):
    '''
    '''
    def __call__(self, atom):
        return OEGetResidueIndex(atom) == OEResidueIndex_HOH or OEAtomGetResidue(atom).GetName() == 'DOD'

class OEIsMetal(OEUnaryAtomPred):
    '''
    '''
    def __call__(self, atom):
        return atom.IsMetal()

class OEIsStdProteinResidue(OEUnaryAtomPred):
    '''
    Returns True is the atom is part of one of the standard amino acids.
    '''
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in AMINOACIDS

class OEIsUnknownResidue(OEUnaryAtomPred):
    '''
    Returns True is the atom is part of one of the standard amino acids.
    '''
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in UNKNOWN_RES

class OEIsNucleotide(OEUnaryAtomPred):
    '''
    Returns True is the atom is part of one of the standard amino acids.
    '''
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in NUCLEOTIDES

class OEIsUnknownElement(OEUnaryAtomPred):
    '''
    '''
    def __call__(self, atom):
        return atom.GetAtomicNum() == 0

OENotIsWater = OENotAtom(OEIsWater())
OEIsInGLY = OEHasResidueIndex(OEResidueIndex_GLY)
OENotBackboneAtom = OENotAtom(OEIsBackboneAtom())

# shortcut to include all polymer atoms
OEIsPolymerAtom = OEAtomBinaryAndIntData(('entity_type_bm', 60))

def file_to_mol(path):
    '''
    '''
    if os.path.exists(path):
        ifs = oemolistream()
        ifs.open(path)

        try: molecule = ifs.GetOEGraphMols().next()
        except StopIteration: return None
        
        return molecule

def pdb_to_mol(path):
    '''
    '''
    if not os.path.exists(path): raise IOError('No such PDB structure: {0}'.format(path))

    # PDB STRUCTURE READER
    ifs = oemolistream()
    ifs.SetFlavor(OEFormat_PDB, OEIFlavor_Generic_Default | OEIFlavor_PDB_ALTLOC | OEIFlavor_PDB_Default | OEIFlavor_PDB_DATA)

    ifs.open(str(path))

    # BY DEFAULT, NMR MODELS ARE TREATED AS SEQUENTIAL MOLECULES, WHICH MEANS
    # ONLY THE FIRST MODEL IS CONSIDERED HERE
    structure = ifs.GetOEGraphMols().next()

    OEAssignAromaticFlags(structure)

    return structure

def smiles_to_mol(smiles, canon=False, strict=True):
    '''
    Returns an OEGraphMol for the given SMILES string.
    '''
    if isinstance(smiles, str):
        molecule = OEGraphMol()
        OEParseSmiles(molecule, smiles, canon, strict)

        return molecule

def standardise_smiles(smiles, isomeric=True, reset_charges=True, kekule=False, neutral=False):
    '''
    Returns the canonical SMILES string of a molecule.
    '''
    smiflag = OESMILESFlag_Canonical

    molecule = OEGraphMol()
    OEParseSmiles(molecule, smiles)

    # RESETS FORMAL CHARGES
    if reset_charges:
        for atom in molecule.GetAtoms(OEIsPolar()):
            atom.SetImplicitHCount(0)
            atom.SetFormalCharge(0)

        OEAssignImplicitHydrogens(molecule)
        OEAssignFormalCharges(molecule)

    OEFindRingAtomsAndBonds(molecule)
    OEAssignAromaticFlags(molecule,OEAroModelOpenEye)

    # APPLIES A FORMAL CHARGE MODEL BASED ON PHYSIOLOGICAL PH
    if neutral: OESetNeutralpHModel(molecule)

    # CREATE SMILES WITH STEREO INFORMATION
    if isomeric:
        smiflag|=OESMILESFlag_ISOMERIC

        OEPerceiveChiral(molecule)

        # CLEAR STEREO AT NONCHIRAL ATOMS
        for atom in molecule.GetAtoms(OEAndAtom(OEHasAtomStereoSpecified(), OENotAtom(OEIsChiralAtom()))):
            nbrs = []

            for bond in atom.GetBonds(): nbrs.append(bond.GetNbr(atom))
            atom.SetStereo(nbrs, OEAtomStereo_Tetrahedral, OEAtomStereo_Undefined)

        # CLEAR STEREO AT NONCHIRAL BONDS
        for bond in molecule.GetBonds(OEAndBond(OEHasOrder(2), OEAndBond(OEHasBondStereoSpecified(), OENotBond(OEIsChiralBond())))):
            bond.SetStereo([bond.GetBgn(),bond.GetEnd()], OEBondStereo_CisTrans, OEBondStereo_Undefined)

    if kekule:
        for bond in molecule.GetBonds(OEIsAromaticBond()): bond.SetIntType(5)

        OECanonicalOrderAtoms(molecule)
        OECanonicalOrderBonds(molecule)
        OEClearAromaticFlags(molecule)
        OEKekulize(molecule)

    ism = OECreateSmiString(molecule,smiflag)

    return ism

def mol_to_smiles(molecule, isomeric=True, reset_charges=True, kekule=False, from3d=False, neutral=False):
    '''
    Returns the canonical SMILES string of a molecule.
    '''
    smiflag = OESMILESFlag_Canonical

    # RESETS FORMAL CHARGES
    if reset_charges:
        for atom in molecule.GetAtoms(OEIsPolar()):
                atom.SetImplicitHCount(0)
                atom.SetFormalCharge(0)

        OEAssignImplicitHydrogens(molecule)
        OEAssignFormalCharges(molecule)

    OEFindRingAtomsAndBonds(molecule)
    OEAssignAromaticFlags(molecule,OEAroModelOpenEye)

    # APPLIES A FORMAL CHARGE MODEL BASED ON PHYSIOLOGICAL PH
    if neutral: OESetNeutralpHModel(molecule)

    # GET STEREO INFORMATION FROM 3D COORDINATES
    if from3d: OE3DToInternalStereo(molecule)

    # CREATE SMILES WITH STEREO INFORMATION
    if isomeric:
        smiflag|=OESMILESFlag_ISOMERIC

        OEPerceiveChiral(molecule)

        # CLEAR STEREO AT NONCHIRAL ATOMS
        for atom in molecule.GetAtoms(OEAndAtom(OEHasAtomStereoSpecified(), OENotAtom(OEIsChiralAtom()))):
            nbrs = []

            for bond in atom.GetBonds(): nbrs.append(bond.GetNbr(atom))
            atom.SetStereo(nbrs, OEAtomStereo_Tetrahedral, OEAtomStereo_Undefined)

        # CLEAR STEREO AT NONCHIRAL BONDS
        for bond in molecule.GetBonds(OEAndBond(OEHasOrder(2), OEAndBond(OEHasBondStereoSpecified(), OENotBond(OEIsChiralBond())))):
            bond.SetStereo([bond.GetBgn(),bond.GetEnd()], OEBondStereo_CisTrans, OEBondStereo_Undefined)

    if kekule:
        for bond in molecule.GetBonds(OEIsAromaticBond()): bond.SetIntType(5)

        OECanonicalOrderAtoms(molecule)
        OECanonicalOrderBonds(molecule)
        OEClearAromaticFlags(molecule)
        OEKekulize(molecule)

    ism = OECreateSmiString(molecule,smiflag)

    return ism

def mol_to_oeb(molecule):
    '''
    '''
    oss = oemolostream()
    oss.SetFormat(OEFormat_OEB)

    oss.openstring()
    OEWriteMolecule(oss,molecule)
    oeb = oss.GetString()

    return oeb

def mol_to_pdb(molecule):
    '''
    '''
    oss = oemolostream()
    oss.SetFormat(OEFormat_PDB)

    oss.openstring()
    OEWriteMolecule(oss,molecule)
    pdb = oss.GetString()

    return pdb

def get_subsetmol(molecule, predicate):
    '''
    Returns the subset of a molecule defined by a given predicate as new molecule.
    '''
    subsetmol = OEMol()

    subsetatoms = molecule.GetAtoms(predicate)
    member = OEIsAtomMember(subsetatoms)

    OESubsetMol(subsetmol, molecule, member, False, False)

    return subsetmol

def reset_charges(molecule):
    '''
    '''
    for atom in molecule.GetAtoms(OEIsPolar()):
        atom.SetImplicitHCount(0)
        atom.SetFormalCharge(0)

    OEAssignImplicitHydrogens(molecule)
    OEAssignFormalCharges(molecule)
    OEAssignAromaticFlags(molecule)

    return molecule

def has_undef_stereo(molecule):
    '''
    '''
    OEPerceiveChiral(molecule)

    predicate = OEAndAtom(OENotAtom(OEHasAtomStereoSpecified()), OEIsChiralAtom())

    return OECount(molecule, predicate) > 0

def get_state_with_largest_avg_occupancy(structure):
    '''
    '''
    altstruct = OEGraphMol()

    # KEEP ONLY THE LOCATION STATE WITH THE LARGEST AVERAGE OCCUPANCY
    altlocfactory = OEAltLocationFactory(structure)
    altlocfactory.MakeCurrentAltMol(altstruct)

    return altstruct

def get_solvent_accessible_surface_areas(molecule):
    '''
    '''
    # ASSIGN VDW RADII / REQUIRED FOR ASA
    OEAssignBondiVdWRadii(molecule)

    # MAKE SOLVENT-ACCESSIBLE SURFACE
    surface = OESurface()
    OEMakeAccessibleSurface(surface, molecule, 0.5, 1.4)

    # IDENTIFY POLAR ATOMS EXPOSED ON THE SURFACE AND CREATE CLIQUE
    for i in xrange(surface.GetNumVertices()):
        atom = molecule.GetAtom(OEHasAtomIdx(surface.GetAtomsElement(i)))

        # WEIRD ERROR WITH CHEMICAL COMPONENT CBW
        if not atom: continue

        # POLAR ATOM IN THIS CONTEXT IS EITHER OXYGEN OR NITROGEN
        if atom.IsOxygen() or atom.IsNitrogen():
            surface.SetVertexCliqueElement(i,1)

    # GET THE TOTAL SURFACE AREA
    asa = OESurfaceArea(surface)

    # GET ONLY THE POLAR SURFACE AREA
    pasa = OESurfaceCliqueArea(surface,1);

    return asa, asa-pasa, pasa