"""
"""

from openeye.oechem import *

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

class OEAtomHasIntData(OEUnaryAtomPred):
    """
    """
    def __init__(self, pair):
        OEUnaryAtomPred.__init__(self)
        self.pair = pair

    def __call__(self, atom):
        return atom.GetIntData(self.pair[0]) == self.pair[1]

    def CreateCopy(self):
        return OEAtomHasIntData(self.pair).__disown__()

class OEAtomHasAnyIntData(OEUnaryAtomPred):
    """
    """
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
    """
    """
    def __init__(self, pair):
        OEUnaryAtomPred.__init__(self)
        self.pair = pair

    def __call__(self, atom):
        return atom.GetIntData(self.pair[0]) & self.pair[1] == 0

    def CreateCopy(self):
        return OEAtomBinaryNotIntData(self.pair).__disown__()

class OEHasResidueIndex(OEUnaryAtomPred):
    """
    """
    def __init__(self, index):
        OEUnaryAtomPred.__init__(self)
        self.index = index

    def __call__(self, atom):
        return OEGetResidueIndex(atom) == self.index

    def CreateCopy(self):
        return OEHasResidueIndex(self.index).__disown__()

class OEIsResidue(OEUnaryAtomPred):
    """
    Returns True if the atom belongs to the given residue.
    """
    def __init__(self, residue):
        OEUnaryAtomPred.__init__(self)
        self.residue = residue

    def __call__(self, atom):
        return OESameResidue(self.residue, OEAtomGetResidue(atom))

    def CreateCopy(self):
        return OEIsResidue(self.residue).__disown__()

class OEIsStdAtom(OEUnaryAtomPred):
    """
    Returns True if the atom is a standard organic element.
    """
    def __call__(self, atom):
        return atom.GetAtomicNum() in STD_ELEMENTS

class OEIsHetatm(OEUnaryAtomPred):
    """
    Returns True if the atom is a standard organic element.
    """
    def __call__(self, atom):
        return OEAtomGetResidue(atom).IsHetAtom()

class OEIsWater(OEUnaryAtomPred):
    """
    """
    def __call__(self, atom):
        return OEGetResidueIndex(atom) == OEResidueIndex_HOH or OEAtomGetResidue(atom).GetName() == 'DOD'

class OEIsMetal(OEUnaryAtomPred):
    """
    """
    def __call__(self, atom):
        return atom.IsMetal()

class OEIsStdProteinResidue(OEUnaryAtomPred):
    """
    Returns True is the atom is part of one of the standard amino acids.
    """
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in AMINOACIDS

class OEIsUnknownResidue(OEUnaryAtomPred):
    """
    Returns True is the atom is part of one of the standard amino acids.
    """
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in UNKNOWN_RES

class OEIsNucleotide(OEUnaryAtomPred):
    """
    Returns True is the atom is part of one of the standard amino acids.
    """
    def __call__(self, atom):
        return OEGetResidueIndex(atom) in NUCLEOTIDES

class OEIsUnknownElement(OEUnaryAtomPred):
    """
    """
    def __call__(self, atom):
        return atom.GetAtomicNum() == 0

OENotIsWater = OENotAtom(OEIsWater())
OEIsInGLY = OEHasResidueIndex(OEResidueIndex_GLY)
OENotBackboneAtom = OENotAtom(OEIsBackboneAtom())

# shortcut to include all polymer atoms
OEIsPolymerAtom = OEAtomBinaryAndIntData(('entity_type_bm', 60))