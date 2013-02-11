from openeye.oechem import (OEGetBondiVdWRadius, OEGetAngle, OEGetDistance,
                            OEGetResidueIndex, OEGetSingleBondNeighbor, OEIsHydrogen,
                            OEResidueIndex_HOH, OEAtomGetResidue)

from credovi import app

contact_types = app.config['contact types']

# a hydrogen bond has some features of covalent bonding, namely interatomic distances shorter than
# the sum of the van der waals radii. since we do not normally know the locations of hydrogens, we
# simply use the following formula as distance cutoff: the distance between the donor (d) and the
# acceptor (a) must be less than or equal to the distance between d and the hydrogen (h) + the vdw
# radius of h and the vdw radius of a. this should be evidence of the partial covalent-like character.
VDW_HYDROGEN            = OEGetBondiVdWRadius(1)
BOND_LENGTH_HYDROGEN    = 1.0

# the bond length of hydrogen and its vdw radius will always be used together
HYDROGEN_CONTRIB        = BOND_LENGTH_HYDROGEN + VDW_HYDROGEN

# cutoff for water-water hydrogen bonds
HBOND_HOH_DIST          = HYDROGEN_CONTRIB + OEGetBondiVdWRadius(8)

# compensation factor in angstrom that is added to contact types calculations to account for xtal uncertainty
VDW_COMP_FACTOR         = 0.1

def is_hbond(structure, atom_bgn, atom_end, distance):
    """
    """
    def _is_hbond(donor, acceptor):
        """
        """
        for hydrogen in donor.GetAtoms(OEIsHydrogen()):
            if OEGetDistance(structure, hydrogen, acceptor) <= VDW_HYDROGEN + acceptor.GetRadius() + VDW_COMP_FACTOR:

                # angle between donor, hydrogen and acceptor
                if OEGetAngle(structure, donor, hydrogen, acceptor) >= contact_types['hbond']['angle rad']:
                    return 1

        return 0

    IS_HBOND = 0

    # one of the atoms is water - pointless to use hydrogen
    if OEGetResidueIndex(atom_bgn) == OEResidueIndex_HOH:
        if (atom_end.GetIntData('hbond donor') or atom_end.GetIntData('hbond acceptor')) and distance <= HBOND_HOH_DIST:
            IS_HBOND = 1

    elif OEGetResidueIndex(atom_end) == OEResidueIndex_HOH:
        if (atom_bgn.GetIntData('hbond donor') or atom_bgn.GetIntData('hbond acceptor')) and distance <= HBOND_HOH_DIST:
            IS_HBOND = 1

    # both atoms are not water
    # we use the OR boolean expression to make sure that a found hydrogen bond
    # is not overwritten by 0 again!
    else:

        # atom bgn is donor
        if atom_bgn.GetIntData('hbond donor') and atom_end.GetIntData('hbond acceptor'):
            IS_HBOND = IS_HBOND or _is_hbond(atom_bgn, atom_end)

        # atom end is donor
        if atom_end.GetIntData('hbond donor') and atom_bgn.GetIntData('hbond acceptor'):
            IS_HBOND = IS_HBOND or _is_hbond(atom_end, atom_bgn)

    return IS_HBOND

def is_weak_hbond(structure, atom_bgn, atom_end, distance):
    """
    """
    def _is_weak_hbond(donor, acceptor):
        """
        """
        for hydrogen in donor.GetAtoms(OEIsHydrogen()):

            # check for bond-like character
            if OEGetDistance(structure, hydrogen, acceptor) <= VDW_HYDROGEN + acceptor.GetRadius() + VDW_COMP_FACTOR:
                if OEGetAngle(structure, donor, hydrogen, acceptor) >= contact_types['weak hbond']['angle rad']:
                    return 1

        return 0

    def _is_halogen_weak_hbond(donor, halogen):
        """
        Halogens can acts as weak hydrogen bond acceptors in a head-on orientation.
        """
        # atom attached to halogen - important for the identification of the side-on orientation
        nbr = OEGetSingleBondNeighbor(halogen)

        for hydrogen in donor.GetAtoms(OEIsHydrogen()):

            # check for bond-like character
            if OEGetDistance(structure, halogen, hydrogen) <= VDW_HYDROGEN + halogen.GetRadius() + VDW_COMP_FACTOR:

                # angle between the neighbor, the halogen acceptor and the hydrogen
                if contact_types['weak hbond']['cx angle min rad'] <= OEGetAngle(structure, nbr, halogen, hydrogen) <= contact_types['weak hbond']['cx angle max rad']:
                    return 1

        return 0

    IS_WEAK_HBOND = 0

    # ATOM IS ACCEPTOR / HETATM IS CARBON
    if (atom_bgn.GetIntData('hbond acceptor') and atom_end.GetIntData('weak hbond donor')):
        IS_WEAK_HBOND = IS_WEAK_HBOND or _is_weak_hbond(atom_end, atom_bgn)

    # HETATM IS ACCEPTOR / ATOM IS CARBON
    if (atom_bgn.GetIntData('weak hbond donor') and atom_end.GetIntData('hbond acceptor')):
        IS_WEAK_HBOND = IS_WEAK_HBOND or _is_weak_hbond(atom_bgn, atom_end)

    # ATOM_BGN IS HALOGEN WEAK ACCEPTOR
    if atom_bgn.GetIntData('weak hbond acceptor') and atom_bgn.IsHalogen() and (atom_end.GetIntData('hbond donor') or atom_end.GetIntData('weak hbond donor')):
        IS_WEAK_HBOND = IS_WEAK_HBOND or _is_halogen_weak_hbond(atom_end, atom_bgn)

    # ATOM_END IS HALOGEN WEAK ACCEPTOR
    if atom_end.GetIntData('weak hbond acceptor') and atom_end.IsHalogen() and (atom_bgn.GetIntData('hbond donor') or atom_bgn.GetIntData('weak hbond donor')):
        IS_WEAK_HBOND = IS_WEAK_HBOND or _is_halogen_weak_hbond(atom_bgn, atom_end)

    return IS_WEAK_HBOND

def is_xbond(structure, atom_bgn, atom_end, distance, sum_vdw_radii):
    """
    Halogens can form electrostatic interactions with Lewis-bases (nucleophiles) in a head-on
    orientation.
    """
    def _is_xbond(donor, acceptor):
        """
        """
        #u
        nbr = OEGetSingleBondNeighbor(donor)
        theta = OEGetAngle(structure, nbr, donor, acceptor)

        # angle for the head-on orientation
        if (theta >= contact_types['xbond']['angle theta 1 rad']):
            return 1

        return 0

    IS_XBOND = 0

    if distance <= sum_vdw_radii + VDW_COMP_FACTOR:
        if atom_bgn.GetIntData('xbond donor') and atom_end.GetIntData('xbond acceptor'):
            IS_XBOND = IS_XBOND or _is_xbond(atom_bgn, atom_end)

        if atom_end.GetIntData('xbond donor') and atom_bgn.GetIntData('xbond acceptor'):
            IS_XBOND = IS_XBOND or _is_xbond(atom_end, atom_bgn)

    return IS_XBOND

def is_ionic(atom_bgn, atom_end, distance):
    """
    """
    IS_IONIC = 0

    if distance <= contact_types['ionic']['distance']:
        if atom_bgn.GetIntData('pos ionisable') and atom_end.GetIntData('neg ionisable') : IS_IONIC = 1
        elif atom_bgn.GetIntData('neg ionisable') and atom_end.GetIntData('pos ionisable'): IS_IONIC = 1

    return IS_IONIC

def is_metal_complex(atom_bgn, atom_end, distance):
    """
    """
    IS_METAL_COMPLEX = 0

    if distance <= contact_types['metal']['distance']:
        if atom_bgn.GetIntData('hbond acceptor') and atom_end.IsMetal(): IS_METAL_COMPLEX = 1
        elif atom_end.GetIntData('hbond acceptor') and atom_bgn.IsMetal(): IS_METAL_COMPLEX = 1

    return IS_METAL_COMPLEX

def is_carbonyl(atom_bgn, atom_end, distance):
    """
    """
    if distance <= contact_types['carbonyl']['distance']:
        if atom_bgn.GetIntData('carbonyl oxygen') and atom_end.GetIntData('carbonyl carbon'): return 1
        elif atom_bgn.GetIntData('carbonyl carbon') and atom_end.GetIntData('carbonyl oxygen'): return 1
    return 0

def is_aromatic(atom_bgn, atom_end, distance):
    """
    """
    if atom_bgn.IsAromatic() and atom_end.IsAromatic() and distance <= contact_types['aromatic']['distance']:
        return 1
    return 0

def is_hydrophobic(atom_bgn, atom_end, distance):
    """
    """
    if atom_bgn.GetIntData('hydrophobe') and atom_end.GetIntData('hydrophobe') and distance <= contact_types['hydrophobic']['distance']:
        return 1
    return 0
