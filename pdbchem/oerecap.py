#!/usr/bin/env python
from itertools import islice
from eyesopen.oechem import *


# KEEPING LINKERS
REACTIONS = {
    'UREA': '[#7;D2,D3:1]!@[C:2](!@=[#8,#16:3])!@[#7;D2,D3:4]>>[#7:1].[C:2](!@=[#8,#16:3]).[#7:4]',
    'AMIDE': '[C;!$(C([#7])[#7]):1](=!@[#8,#16:2])!@[#7;!D1:3]>>[C:1]=[O:2].[#7:3]',
    'AZO': '[NX2:1]=!@[NX2:2]>>[NX2:1].[NX2:2]',
    'ESTERS': '[C;!D1:1](=!@[O,S:2])!@[O,S;!D1;!$([O,S][*;D1]):3]>>[C:1]=[O,S:2].[O,S:3]',
    'AMINES': '[N;!D1;!$(N-C=[#8]):3](-!@[*;!D1:1])-!@[*;!D1:2]>>[*:1].[N:3].[*:2]', # '[N;!D1;!$(N-C=[#8]):3](-!@[*;!D1:1])-!@[*;!D1:2]>>[*:1].[*:2][N:3]'
    'CYCLIC_AMINES': '[#7;R;D3:1]-!@[*;!D1:2]>>[#7:1].[*:2]',
    'AROMN_ALIC': '[n:1]-!@[C:2]>>[n:1].[C:2]',
    'AROM_AROM': '[c,n:1]-!@[c:2]>>[c,n:1].[c:2]',
    'SULPHONAMIDE': '[#7;D2,D3:1]-!@[S:2](=[O:3])=[O:4]>>[#7:1].[S:2](=[O:3])=[O:4]',
    'ETHERS': '[#6;!D1;!$([#6](=!@[#8,#16])):1]-!@[#8,#16:2]-!@[#6;!D1;!$([#6](=!@[#8,#16])):3]>>[#6:1].[#8,#16:2].[#6:3]',
    'OLEFIN': '[C;!D1:1]=!@[C;!D1:2]>>[C:1].[C:2]',
    'DISULPHIDE': '[S;!D1:1]-!@[S;!D1:2]>>[S:1].[S:2]',
    'RCARBON_RCARBON': '[#6;R:1]-!@[#6;R:2]>>[#6;R:1].[#6;R:2]',
    'PHOSPHATES': '[*;!P;!D1:1][OX2;$([OX2][PX4]):2]>>[*:1].[#8:2]',
    'SULFONYL': '[*;!D1:1][#16X4:2]([*;!D1:3])(=[OX1:4])=[OX1:5]>>[*:1].[#16:2]()(=[O:4])=[O:5].[*:3]',
    'RCARBON_ALIC_RCARBON': '[#6;R:1]-!@[#6;H2:2]-!@[#6;R:3]>>[#6;R:1]-!@[#6:2].[#6;R:3]',
    'NITRO_COMPOUND': '[*;!D1;!#8:1][#7;$([NX3](=O)=O),$([NX3+](=O)[O-]):2]>>[*;!#8:1].[#7:2]',
    'SULFINATO': '[*;!D1:1][S;D3;$(S(=O)(=O)),$(S(=O)([OH])):2]>>[*:1].[S:2]',
    'SULFATES': '[*;!S;!D1:1][OX2;$([OX2][S;$([#16X4](=[OX1])(=[OX1])([OX2])([OX1]))]):2]>>[*:1].[O:2]',
    'SULFONATES': '[*;!S;!D1:1][S;$([#16X4](=[OX1])(=[OX1])([OX1])):2]>>[*:1].[#16:2]',
    'SULFAMOYL': '[*:1][SD4;$(S(=O)(=O)([ND1])):2]>>[*:1].[S:2]',
    'ALDEHYDE': '[CX3H1:2](=[O:3])[#6;R:1]>>[#6:1].[Ch0-1:2]#[O+1:3]' # only aldehydes attached to rings
    #'HALOALKANE': '[#6;!D1;!$([#6]([F,Cl,I,Br])([F,Cl,I,Br])([F,Cl,I,Br]));!$(([F,Cl,I,Br])([F,Cl,I,Br])):1][F,Cl,I,Br;D1:2]>>[#6:1].[F,Cl,I,Br:2]'
}

class RecapHierarchyNode(object):
    '''
    '''
    def __init__(self,mol):
        '''
        '''
        OESuppressHydrogens(mol)

        if mol.GetDimension() != 3: from3d=False
        else: from3d=True

        self.mol = mol
        self.smi = mol_to_smiles(mol, reset_charges=False, from3d=from3d, isomeric=True)
        self.atoms = '+'.join(a.GetName().strip() for a in mol.GetAtoms(OENotAtom(OEHasAtomicNum(0))))
        self.children = {}
        self.parents = {}

    def __repr__(self):
        '''
        '''
        return "<Node(%s)>" % self.smi

    def _gacRecurse(self,res,terminalOnly=False):
        '''
        '''
        for smi,child in self.children.iteritems():
            if not terminalOnly or not len(child.children):
                res[smi] = child
            child._gacRecurse(res,terminalOnly=terminalOnly)

    def getAllChildren(self):
        '''
        '''
        " returns a dictionary, keyed by SMILES, of children "
        res = {}
        for smi,child in self.children.iteritems():
            res[smi] = child
            child._gacRecurse(res,terminalOnly=False)
        return res

    def getLeaves(self):
        '''
        returns a dictionary, keyed by SMILES, of leaf (terminal) nodes
        '''
        res = {}
        for smi,child in self.children.iteritems():
            if not len(child.children):
                res[smi] = child
            else:
                child._gacRecurse(res,terminalOnly=True)
        return res

def decompose(mol, valence_correction=True, explicit_hydrogens=False, min_fragment_size=0, allow_butyl=False):
    '''
    '''
    atoms = '+'.join(a.GetName().strip() for a in mol.GetAtoms(OENotAtom(OEHasAtomicNum(0))))
    H = RecapHierarchyNode(mol)
    actives = {(H.smi, atoms):H}
    nodes = {(H.smi, atoms):H}

    if mol.GetDimension() != 3:
        f  rom3d = False
    else:
        from3d = True
    cycles = 0
    BREAK_CYCLE = False
    while actives and not BREAK_CYCLE:
        cycles += 1
        nodeSMI, nodeAtoms = actives.keys()[0]
        node = actives.pop((nodeSMI, nodeAtoms))
        for name, reaction in REACTIONS.items():
            libgen = OELibraryGen(reaction)
            libgen.SetValenceCorrection(valence_correction)
            libgen.SetExplicitHydrogens(explicit_hydrogens)
            libgen.SetStartingMaterial(node.mol, 0)
            
            for product in libgen.GetProducts():
                fragments = []
                VALID_REACTION = True
                # FRAGMENT MOLECULE
                partcount, partlist = OEDetermineComponents(product)
                pred = OEPartPredAtom(partlist)

                for i in range(partcount):
                    fragment = OEGraphMol()
                    pred.SelectPart(i+1)
                    OESubsetMol(fragment, product, pred)
                    num_dummys = OECount(fragment, OEHasAtomicNum(0))
                    num_atoms = OECount(fragment, OENotAtom(OEHasAtomicNum(0)))
                    num_carbons = OECount(fragment, OEIsCarbon())

                    if num_dummys == fragment.NumAtoms():
                        continue
                    # CHECK FOR FORBIDDEN FRAGMENTS AND IF FRAGMENT HAS MINIMUM NUMBER OF HEAVY ATOMS
                    if (not allow_butyl and num_carbons <= 3 and num_atoms == num_carbons) or (num_atoms < min_fragment_size):
                        VALID_REACTION = False
                        break
                    fragments.append(fragment)
                    
                if VALID_REACTION:
                    for fragment in fragments:
                        pSMI = mol_to_smiles(fragment, from3d=from3d, isomeric=True)
                        if pSMI == '*N': 
                            print name, nodeSMI
                        pAtoms = '+'.join(a.GetName().strip() for a in fragment.GetAtoms(OENotAtom(OEHasAtomicNum(0))))

                        if not nodes.has_key((pSMI,pAtoms)):
                            pNode = RecapHierarchyNode(fragment)
                            node.children[(pSMI,pAtoms)] = pNode
                            actives[(pSMI,pAtoms)] = pNode
                            nodes[(pSMI,pAtoms)] = pNode
                        else:
                            pNode = nodes[(pSMI,pAtoms)]
                            node.children[(pSMI,pAtoms)] = pNode

        if cycles > 500: BREAK_CYCLE = True

    # DO ONE PASS CLEAVAGE
    if BREAK_CYCLE:

        atoms = '+'.join(a.GetName().strip() for a in mol.GetAtoms(OENotAtom(OEHasAtomicNum(0))))
        H = RecapHierarchyNode(mol)
        actives = {(H.smi, atoms):H}
        nodes = {(H.smi, atoms):H}

        while actives:
            nodeSMI, nodeAtoms = actives.keys()[0]
            node = actives.pop((nodeSMI, nodeAtoms))

            for name, reaction in REACTIONS.items():
                product = OEGraphMol(node.mol)
                umr = OEUniMolecularRxn(reaction)
                umr(product)

                fragments = []
                VALID_REACTION = True

                # FRAGMENT MOLECULE
                partcount,partlist = OEDetermineComponents(product)
                pred = OEPartPredAtom(partlist)

                # NO REACTION
                if partcount == 1: continue

                for i in range(partcount):
                    fragment = OEGraphMol()
                    pred.SelectPart(i+1)
                    OESubsetMol(fragment, product, pred)

                    num_dummys = OECount(fragment, OEHasAtomicNum(0))
                    num_atoms = OECount(fragment, OENotAtom(OEHasAtomicNum(0)))
                    num_carbons = OECount(fragment, OEIsCarbon())

                    if num_dummys == fragment.NumAtoms(): continue

                    # CHECK FOR FORBIDDEN FRAGMENTS AND IF FRAGMENT HAS MINIMUM NUMBER OF HEAVY ATOMS
                    if (not allow_butyl and num_carbons <= 3 and num_atoms == num_carbons) or (num_atoms < min_fragment_size):
                        VALID_REACTION = False
                        break

                    fragments.append(fragment)

                if VALID_REACTION:
                    for fragment in fragments:
                        pSMI = mol_to_smiles(fragment, from3d=from3d, isomeric=True)
                        pAtoms = '+'.join(a.GetName().strip() for a in fragment.GetAtoms(OENotAtom(OEHasAtomicNum(0))))

                        if not nodes.has_key((pSMI,pAtoms)):
                            pNode = RecapHierarchyNode(fragment)
                            node.children[(pSMI,pAtoms)] = pNode
                            actives[(pSMI,pAtoms)] = pNode
                            nodes[(pSMI,pAtoms)] = pNode

                        else:
                            pNode = nodes[(pSMI,pAtoms)]
                            node.children[(pSMI,pAtoms)] = pNode

    return H