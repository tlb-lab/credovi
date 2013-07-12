import numpy as np
from numpy.linalg import norm
from openeye.oechem import OEAtomGetResidue

def normalized(vector):
    """
    Returns a normalized copy of the given vector.
    """
    return vector / norm(vector)

class AromaticRing(object):
    """
    """
    def __init__(self, number, atoms):

        # an aromatic ring must have a minimum number of atoms
        if len(atoms) < 4:
            raise ValueError("cannot find normal: ring has less than 4 atoms.")

        # get parent structure of the atoms
        structure = atoms[0].GetParent()

        self.residue = OEAtomGetResidue(atoms[0])
        self.number = number
        self.atoms = atoms
        self.pdb = structure.GetTitle()
        self.pdb_chain_id = self.residue.GetChainID()
        self.res_num = self.residue.GetResidueNumber()
        self.res_name = self.residue.GetName().strip()
        self.ins_code = self.residue.GetInsertCode().strip()

        self.name = "/{0}//{1}/{2}`{3}{4}/AR{5}".format(self.pdb, self.pdb_chain_id,
                                                        self.res_name, self.res_num,
                                                        self.ins_code, number)

        # get the coordinates of all the atoms
        coords = np.array([structure.GetCoords(atom) for atom in atoms])

        # centroid is simply the average of coordinates
        self.centroid = coords.mean(0)

        # get two vectors inside the aromatic ring
        AB = coords[0] - coords[1]
        AC = coords[0] - coords[2]

        # the normal is now simply the normalized cross product
        self.normal = normalized(np.cross(AB,AC))

    def __eq__(self, other):
        """
        """
        if isinstance(other, AromaticRing):
            return self.name == other.name
        
        return False

    def __repr__(self):
        """
        """
        return self.name

    def __sub__(self, other):
        """
        Perform vector subtraction.
        """
        return self.centroid - other.centroid

    def distance(self, other):
        """
        Returns the inter-centroid distance between two aromatic rings.
        """
        return norm(self.centroid - other.centroid)

    def angle(self, other, degrees=False, signed=False):
        """
        """
        # other is AromaticRing instance
        if isinstance(other, AromaticRing):
            cosangle = np.dot(self.normal, other.normal) / (norm(self.normal) * norm(other.normal))

        # other is numpy array
        elif isinstance(other, np.ndarray):
            cosangle = np.dot(self.normal, other) / (norm(self.normal) * norm(other))

        # get the angle as radians
        rad = np.arccos(cosangle)

        if not degrees: return rad

        # convert radians into degrees
        else:

            # convert into a signed angle
            if signed: rad = rad -np.pi if rad > np.pi / 2 else rad

            # return degrees
            return rad * 180 / np.pi

if __name__ == '__main__':
    pass

