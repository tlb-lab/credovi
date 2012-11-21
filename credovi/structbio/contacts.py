"""
"""

import os
from math import sqrt
from cStringIO import StringIO
from itertools import combinations, compress, groupby

import tablib
from eyesopen.oechem import *

from credovi import __path__, app
from credovi.structbio import AromaticRing
from credovi.structbio import structure as struct
from credovi.structbio.interactions import (is_hbond, is_weak_hbond, is_xbond,
                                            is_aromatic, is_carbonyl, is_metal_complex,
                                            is_ionic, is_hydrophobic)

def get_contacts(structure, cutoff=5.0, method=OENearestNbrsMethod_Auto):
    """
    Returns an iterator containing the interatomic contact pairs. Only the
    interactions between the given entity and the input structure will be returned
    if an entity is present, otherwise the interactions between all entities
    inside the structure.
    """
    contacts = OEGetNearestNbrs(structure, cutoff, method)

    return contacts

def get_ring_interactions(structure, ringpatterns):
    """
    Returns two tablib data sheets containing information about the aromatic
    rings inside the structure and the interactions that thez form between them.

    :param ringpatterns: List of SMARTS patterns to identify aromatic ring of
                         various sizes
    :returns: Tuple of two tablib data sheets.
    """
    # data set to record aromatic ring details
    ardata = tablib.Dataset(headers=('aromaticring','centroid','normal'),
                            title='aromatic-rings')

    # data set to record ring interaction geometries
    ridata = tablib.Dataset(headers=('bgn','end','distance','dihedral','theta','iota'),
                            title='ring-interactions')

    # get list of all aromatic rings in the structure in the form (hit,[atoms])
    ringatoms = struct.get_aromatic_rings(structure, ringpatterns)

    # get rid of the sequential hit number that is not required here
    _, ringatoms =  zip(*ringatoms)

    # sort aromatic rings before grouping
    # implies that all ring atoms belong to the same residue
    ringatoms = list(ringatoms)
    ringatoms.sort(key=lambda atoms: OEAtomGetResidue(atoms[0]))

    aromaticrings = []
    # group all aromatic rings by residue
    for residue, atomsgroup in groupby(ringatoms,
                                       key=lambda atoms: OEAtomGetResidue(atoms[0])):

        # create a new ring number that is unique in its parent residue
        for i, atoms in enumerate(atomsgroup,1):
            aromaticrings.append(AromaticRing(i, atoms))

    # record the aromatic ring details
    for aromaticring in aromaticrings:
        ardata.append([aromaticring.name,
                       aromaticring.centroid.tolist(),
                       aromaticring.normal.tolist()])

    # get all unique ring interaction pairs
    interactions = combinations(aromaticrings, 2)

    # only keep those within the distance cutoff
    interactions = ((a,b) for a,b in interactions if a.distance(b) <= 6.0)

    # calculate ring interaction geometries
    for a,b in interactions:

        # ignore interactions between aromatic rings of the same residue
        if OESameResidue(a.residue, b.residue): continue

        # calculate the ring interaction geometry
        distance = a.distance(b)
        dihedral = a.angle(b, degrees=True, signed=True)
        theta = a.angle(a-b, True, True)
        iota = b.angle(b-a, True, True)

        ridata.append([a.name, b.name, distance, dihedral, theta, iota])

    return ardata, ridata

def pym(databook, pymfh=None, undefined=True, labels=True):
    """
    Writes a PyMOL script for the given data to visualise interactions.
    """
    # create a string buffer if no file handle was provided
    if not pymfh: pymfh = StringIO()

    # the PyMOL preamble that will be added to the beginning should be in the
    # package data directory
    preamble = os.path.join(__path__[0], 'data', 'preamble.pym')

    # check if preamble actually exist or raise error
    if os.path.isfile(preamble): pymfh.write(open(preamble).read())
    else: raise IOError("cannot open PyMOL preamble: path {} does not exist."
                        .format(preamble))

    # shortcut to the PyMOL dash color configuration
    dashcolors = app.config['pymol']['dashcolor']
    dashradii = app.config['pymol']['dashradius']
    dashgaps = app.config['pymol']['dashgap']
    dashlengths = app.config['pymol']['dashlength']

    for dataset in databook._datasets:

        # create PyMOL command to visualise binary atom interactions
        if dataset.title == 'contacts':

            # template for contact naming convention
            # PyMOL supports regexp-like selection syntax, so this format can
            # easily be used to select a given contact type
            NAME = "{ID1}-{ID2}-{TYPE}-{DIST}"

            # templates for PyMOL commands
            DISTANCE = "distance {NAME}, id {ID1}, id {ID2}\n"
            COLOR = "color {COLOR}, {SELECTION}\n"
            GROUP = "group {GROUP}, {SELECTION}\n"

            # keep track of distance & contact flag combinations
            # important later on to format the dashes and group objects
            found_interactions = set()

            # iterate through our data set and write distance commands for every
            # contact # type
            for row in dataset.dict:

                # row is empty?
                if not row: continue

                # get the atom serials
                bgn, end = row['atom_bgn_serial'], row['atom_end_serial']

                # separate fingerprint into headers and values
                keys, values = zip(*row.items()[5:])

                # compress the list of headers to only set positions on the
                # fingerprint
                interactions = list(compress(keys, values))

                # reformat names of the interaction flag: remove leading 'is_'
                # and any underscores
                interactions = [i[2:].replace('_','') for i in interactions]

                # first item is always the distance flag because one always has
                # to be set
                distflag, contactflags = interactions[0], interactions[1:]

                # undefined contact type, distance probably > 4.5A
                if undefined and not contactflags:
                    name = NAME.format(ID1=bgn, ID2=end, TYPE='undefined', DIST=distflag)
                    distance = DISTANCE.format(NAME=name, ID1=bgn, ID2=end)

                    # write distance command
                    pymfh.write(distance)

                    # add interaction to list of found interactions
                    found_interactions.add(('undefined', distflag))

                # iterate through all identified contact types for this contact
                # multiple interaction types are possible for an interatomic contact,
                # therefore contact types have to be grouped later on to enable the user
                # to switch between them
                else:
                    for contflag in contactflags:
                        name = NAME.format(ID1=bgn, ID2=end, TYPE=contflag, DIST=distflag)
                        distance = DISTANCE.format(NAME=name, ID1=bgn, ID2=end)

                        # write distance command
                        pymfh.write(distance)

                        # add interaction to list of found interactions
                        found_interactions.add((contflag, distflag))

            # color the distances by contact type and format the dashes according to the
            # configuration file
            for contflag, distflag in found_interactions:
                select = "*-*-{}-{}".format(contflag, distflag)

                # write color command to file
                pymfh.write(COLOR.format(COLOR=dashcolors[contflag][distflag],
                                         SELECTION=select))

                # format the dashed lines based on contact and distance flag
                pymfh.write("set {}, {}, {}\n".format('dash_radius', dashradii[contflag][distflag], select))
                pymfh.write("set {}, {}, {}\n".format('dash_gap', dashgaps[contflag][distflag], select))
                pymfh.write("set {}, {}, {}\n".format('dash_length', dashlengths[contflag][distflag], select))

                # hide distance labels if --nolabels option is given
                if not labels: pymfh.write("hide labels, {}\n".format(select))

            # get a unique list of found contact types and group distance objects accordingly
            for contflag in {contflag for contflag, distflag in found_interactions}:

                # group by contact type
                pymfh.write(GROUP.format(GROUP=contflag.upper(),
                                         SELECTION="*-*-{}-*".format(contflag)))

        # write PyMOL commands to visualise aromatic rings
        elif dataset.title == 'aromatic-rings':

            # PyMOL command to create a pseudoatom for the centroid of the ring
            PSEUDOATOM = "pseudoatom {OBJ}, pos={POS}, name={NAME}, chain={CHAIN}, resi={RESNUM}, resn={RESNAME}, color={COLOR}, vdw=0.75\n"

            # iterate through all aromatic rings in this data set
            for row in dataset:

                # row is empty?
                if not row: continue

                ar, centroid, normal = row

                pymfh.write(PSEUDOATOM.format(OBJ=ar.pdb, POS=centroid.tolist(),
                                              NAME="AR" + str(ar.number),
                                              CHAIN=ar.pdb_chain_id, RESNUM=ar.res_num,
                                              RESNAME=ar.res_name,
                                              COLOR=dashcolors['aromatic']['vdw']))

            # show little spheres instead of crosses
            # all centroid (atom) names start with AR
            pymfh.write("hide everything, name AR*\n")
            pymfh.write("show spheres, name AR*\n")

        # write PyMOL commands to visualise ring interactions
        elif dataset.title == 'ring-interactions':

            # iterate through all ring interactions
            for row in dataset:

                # row is empty?
                if not row: continue

                bgn, end = row[:2]

                DISTANCE = "distance {OBJ}, {BGN}, {END}\n"
                pymfh.write(DISTANCE.format(OBJ='RINGINT', BGN=bgn, END=end))

            # color the dashed lines between ring centroids
            pymfh.write("color {}, RINGINT\n".format(dashcolors['aromatic']['vdw']))

            # hide distance labels if --nolabels option is given
            if not labels: pymfh.write("hide labels, RINGINT\n")

    return pymfh

def generate_contacts(structure, **kwargs):
    """
    """
    # these two keyword arguments can be used to limit the contacts to only a
    # specific entity
    pdb_chain_id = kwargs.get('pdb_chain_id', '')
    res_num = kwargs.get('res_num', '')

    # convert the residue number to an integer if valid
    # has to be done separately to avoid bugs when res_num is 0
    if res_num and res_num.isdigit(): res_num = int(res_num)

    # include intramolecular contacts
    intramolecular = kwargs.get('intramolecular', False)

    cutoff = kwargs.get('cutoff', 5.0)
    contact_type_dist_max = kwargs.get('contact_type_dist_max', 4.5)

    # main data book where all data will be stored
    databook = tablib.Databook()

    # create a new and empty data set for this structure
    # headers used for writing contact data
    headers = ('pdb','atom_bgn_serial','atom_end_serial','distance','is_intramolecular',
               'is_clash','is_covalent','is_vdw_clash','is_vdw','is_proximal',
               'is_hbond','is_weak_hbond','is_xbond','is_ionic','is_metal_complex',
               'is_aromatic','is_hydrophobic','is_carbonyl')

    # create a dataset for interatomic contacts
    ctdata = tablib.Dataset(headers=headers, title='contacts')

    # make sure that aromaticity is perceived
    OEAssignAromaticFlags(structure)

    # Determine hyb of all atoms in structure could be useful later on
    OEAssignHybridization(structure)

    # assign Bondi vdw radii
    OEAssignBondiVdWRadii(structure)

    # assign tripos atom names to all atoms
    OETriposAtomTypeNames(structure)

    # set credo atom type flags to all atoms
    struct.set_atom_type_flags(structure)

    # identify all disconnected components including solvent and assign unique
    # entity serial numbers to each, including single atoms/ions.
    numparts, pred = struct.identify_disconnected_components(structure)

    # identify the surface atoms of all entities in asymmetric structure
    # this step has to be done AFTER entity_serials are assigned!
    struct.identify_surface_atoms(structure)

    # label all water molecules as solvents
    for atom in structure.GetAtoms(OEIsWater()):
        atom.SetIntData('is_exposed', 1)

    # calculate ring-interaction geometries
    if kwargs.get('ri'):
        ardata, ridata = get_ring_interactions(structure,
                                               app.config['atom types']['aromatic'].values())
        databook.add_sheet(ardata)
        databook.add_sheet(ridata)

    # get all by all contacts / includes covalently bound neighbour atoms
    contacts = get_contacts(structure, cutoff=app.config['cutoffs']['cutoff'])

    # protonate structure if necessary in order to calculate hbond angles
    # has to be done AFTER determining the contacts to avoid hydrogens showing
    # up in the contact pairs
    if not OEHasExplicitHydrogens(structure):
        OEAddExplicitHydrogens(structure, False, True)
        OESet3DHydrogenGeom(structure)

    # iterate through contacts
    for contact in contacts:
        '''
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
        '''
        # initialize structural interaction fingerprint
        SIFt = [0] * 13
        IS_INTRAMOLECULAR = False

        # get interacting atoms
        if contact.GetBgn() < contact.GetEnd():
            atom_bgn, atom_end = contact.GetBgn(), contact.GetEnd()
        else:
            atom_bgn, atom_end = contact.GetEnd(), contact.GetBgn()

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
            if atom_end.GetIntData('entity_serial') == 0:
                continue

            # second atom belongs to entity / ignore contacts with non-exposed entity atoms
            if atom_end.GetIntData('is_exposed') == 0: continue

        # get the residues
        res_bgn, res_end = OEAtomGetResidue(atom_bgn), OEAtomGetResidue(atom_end)

        # check for selection
        if pdb_chain_id and not res_num:
            if res_bgn.GetChainID() != pdb_chain_id and res_end.GetChainID() != pdb_chain_id:
                continue

        elif pdb_chain_id and res_num:
            if (res_bgn.GetChainID() != pdb_chain_id and res_bgn.GetResidueNumber() != res_num and
                res_end.GetChainID() != pdb_chain_id and res_end.GetResidueNumber() != res_num):
                continue

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
            SIFt[5] = is_hbond(structure, atom_bgn, atom_end, distance)
            SIFt[6] = is_weak_hbond(structure, atom_bgn, atom_end, distance)
            SIFt[7] = is_xbond(structure, atom_bgn, atom_end, distance, sum_vdw_radii)
            SIFt[8] = is_ionic(atom_bgn, atom_end, distance)
            SIFt[9] = is_metal_complex(atom_bgn, atom_end, distance)
            SIFt[10] = is_aromatic(atom_bgn, atom_end, distance)
            SIFt[11] = is_hydrophobic(atom_bgn, atom_end, distance)
            SIFt[12] = is_carbonyl(atom_bgn, atom_end, distance)

        ### write contact details to file / atom serial is used to identify atom

        row = ['\N' for i in range(5)]

        row[0] = structure.GetTitle() # PDB
        row[1] = res_bgn.GetSerialNumber()
        row[2] = res_end.GetSerialNumber()
        row[3] = distance

        # write boolean as 0 or 1
        row[4] = 1 if IS_INTRAMOLECULAR else 0

        # add SIFt to contact info
        row.extend(SIFt)

        # add this row to our data set for this structure
        ctdata.append(row)

    # add contact sheet to our main data book
    databook.add_sheet(ctdata)

    return databook
