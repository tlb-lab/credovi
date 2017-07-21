#!/usr/bin/env python
import numpy as np
from itertools import count, product, permutations, combinations
from collections import deque


PHARMACOPHORES = ('AC','AP','AR','DO','NE','PO')
RANGES = {'A': (0.0, 4.8),'B': (4.8, 7.2),'C': (7.2, 9.5),'D': (9.5, 11.9),'E': (11.9, 14.3)}
D = sorted([ (k, v[1]**2) for k,v in RANGES.iteritems() ], key=lambda x:x[1])

# MODIFIED BASED ON BLOSUM SUBSTITUTION SCORES
FEATURES  = {
            'MLY':('AC','AP','PO'),
            'HIS':('AC','AR','DO'),    # AC AR DO
            'ASN':('AC','DO'),         # AC DO
            'GLN':('AC','DO'),         # AC DO
            'SER':('AC','DO'),         # AC DO
            'THR':('AC','DO'),         # AP AC DO
            'ASP':('AC','NE'),         # AC NE
            'GLU':('AC','NE'),         # AC NE
            'CGU':('AC','NE'),
            'GLY':('AP',),
            'ALA':('AP',),             # AP
            'ILE':('AP',),             # AP
            'LEU':('AP',),             # AP
            'MET':('AP',),             # AP
            'MSE':('AP',),
            'PRO':('AP',),             # AP
            'VAL':('AP',),             # AP
            'PHE':('AP','AR'),         # AP AR
            'TRP':('AP','AR','DO'),    # AR DO
            'TYR':('AP','AR','DO'),    # AC AR DO
            'CYS':('AP','DO'),         # AP
            'CSO':('AP','DO'),
            'HYP':('AP','DO'),
            'ARG':('DO','PO'),         # AP DO PO
            'LYS':('DO','PO'),         # AP DO PO
            }

REPRESENTATIVES = {'ALA':'CB','ARG':'CZ','ASN':'CG','ASP':'CG','CYS':'SG',
                   'GLN':'CD','GLU':'CD','GLY':'CA','HIS':'NE2','ILE':'CG1',
                   'LEU':'CG','LYS':'CE','MET':'SD','PHE':'CZ','PRO':'CG',
                   'SER':'OG','THR':'CB','TRP':'CD2','TYR':'CZ','VAL':'CB'}

class Devnull(file):
    def write(self, *_): pass
    def __nonzero__(self): return False

class FuzCavTriangle(object):
    """
    """
    def __init__(self,p1,f1,p2,f2,p3,f3):
        """
        """
        self.A = p1
        self.B = p2
        self.C = p3
        self.a = get_range((p2[0]-p3[0])**2+(p2[1]-p3[1])**2+(p2[2]-p3[2])**2)
        self.b = get_range((p1[0]-p3[0])**2+(p1[1]-p3[1])**2+(p1[2]-p3[2])**2)
        self.c = get_range((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
        self.F1 = f1
        self.F2 = f2
        self.F3 = f3

    def isValid(self):
        """
        """
        # CHECK IF ALL TRIPLETS ARE WITHIN THE MAXIMUM ALLOWED DISTANCE (14.3)
        if not all([self.a,self.b,self.c]):
            return False

        # CHECK TRIANGLE INEQUALITY RULE
        if (RANGES[self.a][0] > RANGES[self.b][1] + RANGES[self.c][1]) \
            or (RANGES[self.b][0] > RANGES[self.a][1] + RANGES[self.c][1]) \
            or (RANGES[self.c][0] > RANGES[self.a][1] + RANGES[self.b][1]):
            return False

        return True

    def getFeatures(self):
        """
        """
        return [tuple(sorted(((f1,self.a),(f2,self.b),(f3,self.c)))) for f1 in self.F1 for f2 in self.F2 for f3 in self.F3]

def get_range(distance):
    """
    """
    for rng, cutoff in D:
        if distance <= cutoff:
            return rng

def get_tracker():
    """
    """
    # GENERATE ALL PHARMACOPHORE TRIPLETS / 216
    pharmacophores = list(product(PHARMACOPHORES, repeat=3))

    # GENERATE ALL DISTANCE TRIPLET BINS / 125
    distbins = list(product(sorted(RANGES.keys()), repeat=3))

    # REMOVE DISTANCE TRIPLETS VIOLATING TRIANGLE INEQUALITY RULES / 113
    distbins = [bin for bin in distbins if not (RANGES[bin[0]][0] > RANGES[bin[1]][1] + RANGES[bin[2]][1]
                                                or RANGES[bin[2]][0] > RANGES[bin[0]][1] + RANGES[bin[1]][1]
                                                or RANGES[bin[1]][0] > RANGES[bin[2]][1] + RANGES[bin[0]][1])]

    # CREATE ALL POSSIBLE TRIPLET + TRIPLET COMBINATIONS / 24408
    x = set([tuple(zip(f,d)) for f in pharmacophores for d in distbins])
    triplets = set()

    # ROTATE TRIANGLES, REMOVE ALL REDUNDANT ONCES / 4492
    while x:
        triplet = x.pop()
        d = deque(triplet)

        d.rotate()
        x.discard(tuple(d))

        d.rotate()
        x.discard(tuple(d))

        triplet = tuple(sorted(triplet))
        triplets.add(triplet)

    tracker = dict(zip(triplets, count()))

    return tracker

#def get_fingerprint(features, tracker, bits):
def make_fp(features, tracker):
    """
    """
    #logging.debug("Vector size: %d", len(tracker))
    fingerprint = np.zeros(len(tracker),dtype=np.uint8)
    triplets = combinations(features, 3)

    #print '\n'.join([ str(t) for t in triplets])
    for ((A,F1),(B,F2),(C,F3)) in triplets:
        triangle = FuzCavTriangle(A,F1,B,F2,C,F3) #(A.Coords,F1,B.Coords,F2,C.Coords,F3)

        if not triangle.isValid():
            continue
        for feature in triangle.getFeatures():
            try:
                fingerprint[tracker[feature]] += 1
            except IndexError:
                logging.error("Index %d out of bounds [%d]", tracker[feature], len(fingerprint))

    return fingerprint

def parse_pdbfile(pdbfile, mask_mode=False, **filters):
    atoms     = []
    coords    = []
    pocklines = []
    hetcrds   = []
    pock_mask = []
    het1      = None
    for l,line in enumerate(pdbfile):
        try:
            if line.startswith('ATOM'):
                atomnum = int(line[6:11])
                resnum  = int(line[22:26])
                chain   = line[21]
                atomid  = line[12:16].strip()
                resid3  = line[17:20]
                #resid1  = aa321.get(resid3,'X')  # For later
                if resid3 not in REPRESENTATIVES:
                    logging.warning("Unknown residue found on line %d of file %s: %s-%d", l, pdbfile.name, resid, resnum)
                    if sys.flags.debug:
                        raise KeyError("Unknown residue '%s'" % resid3)
                    else:
                        continue

                if filters.get('residues') and resnum not in filters['residues']:
                    if mask_mode or filters.get('expand'):
                        pock_mask.append(False)
                    else:
                        #logging.debug("Residue filter found. Skipping residue %d (%d).", resnum, atomnum)
                        continue
                elif filters.get('atoms') and atomnum not in filters['atoms']:
                    if mask_mode or filters.get('expand'):
                        pock_mask.append(False)
                    else:
                        #logging.debug("Atom filter found. Atom %d skipped.", atomnum)
                        continue
                elif filters.get('chains') and chain not in filters['chains']:
                    if mask_mode or filters.get('expand'):
                        pock_mask.append(False)
                    else:
                        #logging.debug("Chain filter found. Skipping chain %d (%d)", chain, atomnum)
                        continue
                else:
                    pock_mask.append(True)

                atoms.append((resid3, atomid))
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                if pock_mask[-1] or filters.get('expand'):
                    pocklines.append(line)
            elif line.startswith('HETATM') and filters.get('ligand'):
                resid3  = line[17:20]
                resnum  = line[22:26].strip()
                hetcode = resid3+resnum
                if filters['ligand']=='all' or hetcode.startswith(filters['ligand']) or (filters['ligand']=='first' and hetcode==het1):
                    hetcrds.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    het1 = hetcode
        except (ValueError, IndexError), err:
            logging.error("Error parsing line %d of file %s:\n%s", l, pdbfile.name, line, exc_info=True)
            if sys.flags.debug:
                raise
    #if filters.get('ligand') and len(hetcrds) > 0:
    #    atcoords  = np.array(coords, dtype=float)
    #    hetcoords = np.array(hetcrds, dtype=float)
    #
    #    prox_mask = calc_proximity_mask(atcoords, hetcoords, filters.get('distance',5.0)).astype(np.bool)
    #    logging.info("Atoms within %2.2fA of %s: %d/%d", filters.get('distance',5.0), het1, prox_mask.sum(), len(atcoords))
    #    if mask_mode:
    #        return np.array(atnames, dtype=int), atcoords, np.array(pocklines), prox_mask
    #    else:
    #        return np.array(atnames, dtype=int)[prox_mask], atcoords[prox_mask], np.array(pocklines)[prox_mask], np.array(pock_mask, dtype=bool)[prox_mask]
    #elif filters.get('expand'):
    #    atcoords   = np.array(coords, dtype=float)
    #    pock_mask  = np.array(pock_mask, dtype=bool)
    #
    #    prox_mask  = calc_proximity_mask(atcoords, atcoords[pock_mask], filters['expand']).astype(bool)
    #    logging.info("Pre-defined pocket expanded by %2.2fA. Now includes %d atoms (from %d)", filters['expand'], prox_mask.sum(), pock_mask.sum())
    #    if not mask_mode:
    #        return np.array(atnames, dtype=int)[prox_mask], atcoords[prox_mask], np.array(pocklines)[prox_mask], prox_mask[prox_mask]
    #    else:
    #        pock_mask = prox_mask

    if mask_mode and not filters:
        logging.warning("Using pocket mask mode without any filters does not make sense.")
    return atoms, np.array(coords, dtype=float), pocklines, np.array(pock_mask, dtype=np.bool)


if __name__ != '__main__':
    tracker = get_tracker()
else:
    import sys
    #from progressbar import ProgressBar, Percentage, Bar, SimpleProgress
    import argparse
    import logging
    from os.path import *

    argparser = argparse.ArgumentParser()
    argparser.add_argument('filename', nargs='?', metavar="PDBFILE", type=str, help="Filename of PDB to process", default=None)
    argparser.add_argument('-o', '--out', nargs='?', dest="outfile", type=str, help="File to output signatures into. Default: STDOUT", const='*', default='')
    argparser.add_argument('-p', '--pocket-in',  dest="pocketin", type=argparse.FileType('r'), help="PDB file defining the pocket (mostly useful for pocket mask mode")
    argparser.add_argument('-P', '--pocket-out',  dest="pocketout", type=argparse.FileType('w'), help="File to output pocket to.")
    argparser.add_argument('-M', '--mask', dest='mask_mode', action='store_true', help="Use pocket mask mode, (pocket contacts calculated against all atoms in protein)", default=False)
    argparser.add_argument('-x', '--expand', dest='expand', type=float, help="Expand manually defined pocket (i.e. not ligand based) by X Angstroms", default=0.0)
    argparser.add_argument('-b', '--batch', nargs='?', dest="batch", metavar="BATCHFILE", type=argparse.FileType('r'), help="Batch file with list of files to process. Default: STDIN", const=sys.stdin)
    #argparser.add_argument('-c', '--cutoff-max', dest="cutoff_max", type=float, help="Maximum cutoff limit (Default: 20)",  default=20.0)
    #argparser.add_argument('-s', '--cutoff-step', dest="cutoff_step", type=float, help="Step size (default: 0.2)", default=0.2)
    argparser.add_argument('-t', '--rep-atom', dest="rep_atom", type=str, choices=['ca','dyn'], help="Representative atom to use for signature: all, ca[lpha], dyn[amic]. Default: all", default='all')
    argparser.add_argument('-r', '--residues', nargs='*', dest="residues", metavar="POCKET_RESIDUES", type=str, help="List of residue numbers defining the pocket. Specify ranges as x-y.",  default=[])
    argparser.add_argument('-a', '--atoms', nargs='*', dest="atoms", metavar="POCKET_ATOMS", type=str, help="List of atom numbers defining the pocket. Specify ranges as x-y.",  default=[])
    argparser.add_argument('-l', '--ligand', nargs='?', dest="ligand", type=str, help="Ligand to generate pocket from. Format: LIG123. No argument: First ligand", const='first', default=None)
    argparser.add_argument('-d', '--distance', dest="distance", type=float, help="Threshold distance for ligand pocket generation. Default: 5.0", default=5.0)
    argparser.add_argument('--profile',  dest="profile", action='store_true', help="Run profiling", default=False)
    argparser.add_argument('-L', '--loglevel',dest="loglevel", action='store', choices=['debug','info','warn','error'], help="Logging level to display (default: warnings and errors)", default='warn')

    opt = argparser.parse_args()
    optdict = vars(opt)
    logging.basicConfig(level=getattr(logging, opt.loglevel.upper()))


    pocket = {}
    for elem in ['atoms','residues','ligand','dist','pocketin','expand']:
        if not optdict.get(elem):
            pocket[elem] = None
            continue
        elif elem == 'pocketin':
            pocket['atoms'] = [ int(line[6:11]) for line in opt.pocketin if line.startswith('ATOM') ]
            continue
        elif elem in ['ligand','distance','expand']:
            pocket[elem] = optdict[elem]
            logging.debug(pocket)
            continue
        else:
            pocket[elem] = []
        for e in optdict[elem]:
            if '-' in e:
                first, last = e.split('-')
                pocket[elem].extend([ e for e in range(int(first), int(last)+1) ])
            elif ',' in e:
                pocket[elem].extend(e.split(','))
            else:
                pocket[elem].append(int(e))
        logging.debug("Pocket %s: %s", elem, pocket[elem])

    file_iter = [opt.filename] if opt.filename else opt.batch
    if not file_iter:
        argparser.error("You must provide either a PDB file or enable batch mode (-b)")

    if opt.outfile == '*':
        logging.debug("Outputting signatures of both reps to fuzcav_ca.csv and fuzcav_dyn.csv")
        out_ca  = open('fuzcav_ca.csv', 'a') if opt.rep_atom in ['all','ca'] else Devnull()
        out_rep = open('fuzcav_dyn.csv','a') if opt.rep_atom in ['all','dyn'] else Devnull()
    elif opt.outfile:
        if opt.rep_atom == 'all':
            raw_input("WARNING: You selected to output signature for both representative atoms, but also specified an output file. Both signatures will go into the same place! Press Control-C to cancel, ENTER to continue.")
            out_ca  = out_rep = open(opt.outfile,'w') if opt.rep_atom in ['all','ca']  else Devnull()
        else:
            out_ca  = open(opt.outfile,'w') if opt.rep_atom == 'ca'  else Devnull()
            out_rep = open(opt.outfile,'w') if opt.rep_atom == 'dyn' else Devnull()
    else:
        out_ca  = sys.stdout if opt.rep_atom in ['all','ca']  else Devnull()
        out_rep = sys.stdout if opt.rep_atom in ['all','dyn'] else Devnull()

    tracker = get_tracker()
    for pdb in file_iter:
        pdb = pdb.strip()
        logging.info("Processing %s", pdb)
        try:
            (pdb, pockpdb)  = pdb.split() if len(pdb.split()) > 1 else (pdb, None) # In batch mode, a second entry per line indicates PDB defining pocket, for mask mode or expansion
            pocket['atoms'] = [int(line[6:11]) for line in open(pockpdb) if line.startswith('ATOM')] if (pockpdb and isfile(pockpdb)) else pocket.get('atoms',[])
            with open(pdb) as pdbfile:
                (atoms, coords, pocket_lines, pocket_mask) = parse_pdbfile(pdbfile, mask_mode=opt.mask_mode, **pocket)
            assert len(atoms)==len(pocket_mask)
            assert opt.mask_mode or np.all(pocket_mask)
        except Exception, err:
            logging.error("Error processing file '%s'", pdb, exc_info=True)
            continue
        else:
            logging.info("File %s: %d atoms included. In pocket: %s", pdb, len(atoms), pocket_mask.sum())
            if pocket and opt.pocketout:
                print >> opt.pocketout, ''.join(pocket_lines)
                logging.info("Pocket PDB file %s written. %d lines.", opt.pocketout.name, len(pocket_lines))
        if len(atoms) > 1000 and not opt.mask_mode:
            raw_input("Number of atoms is very large and not running on mask mode. Control-C to cancel, ENTER to continue.")
        elif len(atoms) < 14:
            logging.error("Too few atoms: %s", len(atoms))
            sys.exit(1)

        ca_mask   = np.array([at=='CA' for (res,at) in atoms ], dtype=np.bool)
        rep_mask  = np.array([REPRESENTATIVES.get(res)==at for (res,at) in atoms ], dtype=np.bool)
        features  = np.array([ (coords[i], FEATURES.get(res)) for i,(res,at) in enumerate(atoms) ])

        if opt.profile:
            import pstats, cProfile
            cProfile.runctx("make_fp(features[ca_mask], tracker); make_fp(features[rep_mask], tracker)", globals(), locals(), "fuzcav.prof")
            s = pstats.Stats("fuzcav.prof")
            s.strip_dirs().sort_stats("time").print_stats()
        else:
            logging.debug("Outputting signatures to individual files.")
            root, base = split(pdb)
            base = splitext(base)[0]
            if out_ca:
                fuzfp_ca  = make_fp(features[ca_mask], tracker)
                print >> out_ca, base+',', ','.join([ str(e) for e in fuzfp_ca  ])
                with open(join(root,base+'_fuz_ca.csv'), 'w') as tmpout:
                    print >> tmpout, ','.join([ str(e) for e in fuzfp_ca  ])
            if out_rep:
                fuzfp_rep = make_fp(features[rep_mask], tracker)
                print >> out_rep, base+',', ','.join([ str(e) for e in fuzfp_rep ])
                with open(join(root,base+'_fuz_dyn.csv'), 'w') as tmpout:
                    print >> tmpout, ','.join([ str(e) for e in fuzfp_rep  ])

        #calphas = ((np.array(atom.coords, dtype=float), (fuzcav.FEATURES[res_name]))
        #           for res_name, atom in atoms if atom.atom_name=='CA')
        #
        ## get the representative atom and its features for each residue
        #representatives = ((np.array(atom.coords, dtype=float), (fuzcav.FEATURES[res_name]))
        #                   for res_name, atom in atoms
        #                   if atom.atom_name==fuzcav.REPRESENTATIVES[res_name])



    #if args.progressbar:
    #    bar = ProgressBar(widgets=['Binding Sites: ', SimpleProgress(), ' ',
    #                               Percentage(), Bar()], maxval=len(ligand_ids)).start()
