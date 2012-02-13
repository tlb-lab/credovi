"""
This module is used to generate CREDO data on structures "as they come", i.e.
without using any other external information such as databases.
"""

import os
import sys
from math import sqrt
from itertools import combinations, compress, groupby

import tablib

from credovi import __path__, app
from credovi.structbio import AromaticRing
from credovi.structbio import structure as struct
from credovi.structbio.interactions import (is_hbond, is_weak_hbond, is_xbond,
                                             is_aromatic, is_carbonyl, is_metal_complex,
                                             is_ionic, is_hydrophobic)
from credovi.util.timer import Timer
from credovi.lib.openeye import *

def get_contacts(args, structure):
    """
    """
    timer.start()
            
    ### get all by all contacts / includes covalently bound neighbour atoms
    
    contacts = OEGetNearestNbrs(structure, app.config['cutoffs']['cutoff'],
                                OENearestNbrsMethod_Auto)

    # debug how much time it took to get all contacts
    app.log.debug("all contacts identified in {0:.2f} seconds.".format(timer.elapsed()))

def get_ring_interactions(structure):
    """
    Only used in the offline version.
    """
    # data set to record aromatic ring details
    ardata = tablib.Dataset([],
                            headers=('aromaticring','centroid','normal'),
                            title='aromatic-rings')
    
    # data set to record ring interaction geometries
    ridata = tablib.Dataset([],
                            headers=('bgn','end','distance','dihedral','theta','iota'),
                            title='ring-interactions')
    
    # get list of all aromatic rings in the structure in the form (hit,[atoms])
    ringatoms = struct.get_aromatic_rings(structure,
                                              app.config['atom types']['aromatic'].values())
    
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
        ardata.append([aromaticring, aromaticring.centroid, aromaticring.normal])
    
    # get all unique ring interaction pairs
    interactions = combinations(aromaticrings, 2)    
    
    # only keep those within the distance cutoff
    interactions = ((a,b) for a,b in interactions if a.distance(b) <= 6.0)
    
    # calculate ring interaction geometries
    for a,b in interactions:
        
        # ignore interactions between aromatic rings of the same residue
        if OESameResidue(a.residue, b.residue): continue        
        
        distance = a.distance(b)
        dihedral = a.angle(b, degrees=True, signed=True)
        theta = a.angle(a-b, True, True)
        iota = b.angle(b-a, True, True)
        
        ridata.append([a.name, b.name, distance, dihedral, theta, iota])
            
    return ardata, ridata

def pym(args, databook, filename):
    """
    Writes a PyMOL script for the given data to visualise interactions.
    """
    # if an output directory is given, the PyMOL script will be written
    # to this directory with the same filename but the extension replaced
    # by .pym
    if args.output_dir:
        path = os.path.join(args.output_dir, '{0}.pym'.format(filename))
        pymfh = open(path,'w')     
    
    # if no output directory is given, write to stout
    else: pymfh = sys.stdout    
    
    # copy PyMOL preamble
    preamble = os.path.join(__path__[0], 'data', 'preamble.pym')
    
    if os.path.isfile(preamble): pymfh.write(open(preamble).read())
    else: app.log.warn("unable to load preamble.pym from data directory.")  
    
    # shortcut to the PyMOL dash color configuration
    dashcolors = app.config['pymol']['dashcolor']    
    dashradii = app.config['pymol']['dashradius']
    dashgaps = app.config['pymol']['dashgap']
    dashlengths = app.config['pymol']['dashlength']
    
    for dataset in databook._datasets:

        # create PyMOL command to visualise binary atom interactions
        if dataset.title == 'contacts':
            
            # template for contact naming convention
            # PyMOL supports regexp-like selection syntax, so this format can easily be
            # used to select a given contact type
            NAME = "{ID1}-{ID2}-{TYPE}-{DIST}"
            
            # templates for PyMOL commands
            DISTANCE = "distance {NAME}, id {ID1}, id {ID2}\n"
            COLOR = "color {COLOR}, {SELECTION}\n"     
            GROUP = "group {GROUP}, {SELECTION}\n"    
        
            # keep track of distance & contact flag combinations
            # important later on to format the dashes and group objects    
            found_interactions = set()        
            
            # iterate through our data set and write distance commands for every contact
            # type
            for row in dataset.dict:
                
                # row is empty?
                if not row: continue
                
                # get the atom serials
                bgn, end = row['atom_bgn_serial'], row['atom_end_serial']
        
                # separate fingerprint into headers and values
                keys, values = zip(*row.items()[5:])
                
                # compress the list of headers to only set positions on the fingerprint
                interactions = list(compress(keys, values))
                
                # reformat names of the interaction flag: remove leading 'is_' and any underscores
                interactions = [i[2:].replace('_','') for i in interactions]   
                
                # first item is always the distance flag because one always has to be set
                distflag, contactflags = interactions[0], interactions[1:]       
         
                # undefined contact type, distance probably > 4.5A
                if args.undefined and not contactflags:
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
                if args.nolabels: pymfh.write("hide labels, {}\n".format(select))    
            
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
            if args.nolabels: pymfh.write("hide labels, RINGINT\n")
    
    # close file handle
    pymfh.close()

def do(controller):
    """
    """
    # get the controller command
    cmd = controller.command
    
    # get the command line arguments and options
    args = controller.pargs    
    
    # check if the output directory exists
    if args.output_dir and not os.path.exists(args.output_dir): pass
    
    # OEChem input flavour to use for reading PDB entries
    # this flavour will keep alternate location groups
    OE_PDB_INPUT_FLAVOR = OEIFlavor_Generic_Default | OEIFlavor_PDB_ALTLOC | OEIFlavor_PDB_Default | OEIFlavor_PDB_DATA
        
    # PDB structure reader for gzipped structures
    ifs = oemolistream()
    ifs.SetFlavor(OEFormat_PDB, OE_PDB_INPUT_FLAVOR)
    ifs.Setgz(args.gz)        
    
    # timer to clock functions and parts of the program
    timer = Timer()
    timer.start(controller.Meta.label)       
    
    # iterate through all the files given on the command line
    for path in args.structures.split(','):
        directory, filename = os.path.split(path)
        filename, extension = os.path.splitext(filename)     
        
        # check if the file exists
        if not os.path.isfile(path):
            app.log.error("cannot read structure: file {0} does not exist.".format(path))
            
            # skip to next file
            continue  
    
        # by default, NMR models are treated as sequential molecules, which means
        # only the first model is considered here
        try:
            ifs.open(str(path)) # must not be Unicode
            structure = ifs.GetOEGraphMols().next()
            
            # filename will be used as structure title
            structure.SetTitle(str(filename))
        except StopIteration:
            app.log.error("cannot process structure: unable to parse {0}.".format(path))
            continue
            
        # debug number of atoms in structure
        app.log.debug("structure loaded with {0} atoms.".format(structure.NumAtoms()))

        # main data book where all data will be stored
        databook = tablib.Databook()        
        
        # create a new and empty data set for this structure
        # headers used for writing contact data
        headers = ('pdb','atom_bgn_serial','atom_end_serial','distance','is_intramolecular',
                   'is_clash','is_covalent','is_vdw_clash','is_vdw','is_proximal',
                   'is_hbond','is_weak_hbond','is_xbond','is_ionic','is_metal_complex',
                   'is_aromatic','is_hydrophobic','is_carbonyl')

        ctdata = tablib.Dataset([], headers=headers, title='contacts')
            
        OEAssignAromaticFlags(structure)

        # Determine hyb of all atoms in structure could be useful later on
        OEAssignHybridization(structure)

        # assign Bondi vdw radii
        OEAssignBondiVdWRadii(structure)

        # assign tripos atom names to all atoms
        OETriposAtomTypeNames(structure)        

        # identify all disconnected components including solvent and assign unique
        # entity serial number to each // the function returns the number of
        # identified components as well as an atom idx mapping that can be used
        # to partition the molecule
        numparts, pred = struct.identify_disconnected_components(structure)          
        
        # set credo atom type flags to all atoms
        struct.set_atom_type_flags(structure)        
        
        # identify the surface atoms of all entities in asymmetric structure
        struct.identify_surface_atoms(structure)
        
        # label all water molecules as being exposed
        if args.water:
            for atom in structure.GetAtoms(OEIsWater()):
                atom.SetIntData('is_exposed',1)
                atom.SetIntData('is_solvent',1)        
        
        # calculate ring-interaction geometries
        if args.ri:
            ardata, ridata = get_ring_interactions(structure)
            databook.add_sheet(ardata)
            databook.add_sheet(ridata)

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
        
            IS_COVALENT = False
            IS_INTRAMOLECULAR = False
        
            # get interacting atoms
            if contact.GetBgn() < contact.GetEnd():
                atom_bgn, atom_end = contact.GetBgn(), contact.GetEnd()
            else:
                atom_bgn, atom_end = contact.GetEnd(), contact.GetBgn()
        
            # ignore buried atoms
            if not atom_bgn.GetIntData('is_exposed') or not atom_end.GetIntData('is_exposed'):
                continue
            
            # ignore solvent-solvent contacts
            if atom_bgn.GetIntData('is_solvent') and atom_end.GetIntData('is_solvent'):
                continue
            
            # set a flag for intra-entity contacts
            if atom_bgn.GetIntData('entity_serial') == atom_end.GetIntData('entity_serial'):
                IS_INTRAMOLECULAR = True
        
                # skip if only intermolecular interactions should be recorded
                if not args.intramolecular: continue
        
            # get the residues
            res_bgn, res_end = OEAtomGetResidue(atom_bgn), OEAtomGetResidue(atom_end)
            
            # ignore all intra-residue contacts
            if OESameResidue(res_bgn, res_end): continue
            
            # get interatomic distance
            distance = sqrt(contact.GetDist2())
        
            # ignore intra-molecular contacts that are only separated by three bonds
            if IS_INTRAMOLECULAR and OEGetPathLength(atom_bgn, atom_end, 2): continue

            # get the sum of van der waals radii of both atoms
            sum_vdw_radii = atom_bgn.GetRadius() + atom_end.GetRadius()

            # get the sum of covalent radii
            sum_cov_radii = OEGetCovalentRadius(atom_bgn.GetAtomicNum()) + OEGetCovalentRadius(atom_end.GetAtomicNum()) 

            # use the connection table to identify covalent bonds
            # can identify covalent bonds of unusual length as well
            if atom_end in atom_bgn.GetAtoms():
                SIFt[1] = 1
                IS_COVALENT = True

            # check for atomic clash
            elif distance < sum_cov_radii - 0.1: SIFt[0] = 1

            # check vdw radii
            elif distance < sum_vdw_radii: SIFt[2] = 1
            elif distance <= sum_vdw_radii + 0.1: SIFt[3] = 1 # + vdw comp factor

            # label as proximal
            else: SIFt[4] = 1
        
            # skip this step if atoms are covalently bonded or above the maximum
            # feature contact type distance (4.5) to save time
            if not IS_COVALENT and distance < app.config['cutoffs']['contact type dist max']:

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

            row[0] = '3E5A' # PDB
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

        app.log.debug("all contacts processed in {0:.2f} seconds."
                      .format(timer.elapsed()))
        
        # write contact data to file if necessary
        if args.contacts:
            if args.output_dir:
                path = os.path.join(args.output_dir, '{name}.{ext}'
                                    .format(name=filename, ext=args.format))
                
                # write the data in the specified format
                with open(path,'wb') as fh: fh.write(getattr(ctdata, args.format))
        
        if args.pymol: pym(args, databook, filename)