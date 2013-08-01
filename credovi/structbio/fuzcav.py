#!/usr/bin/env python
from itertools import count, product, permutations, combinations
from collections import deque

from numpy import zeros

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

def make_fp(features, tracker):
    """
    """
    fingerprint = zeros(len(tracker),dtype=int)
    triplets = combinations(features, 3)

    for ((A,F1),(B,F2),(C,F3)) in triplets:
        triangle = FuzCavTriangle(A,F1,B,F2,C,F3)

        if not triangle.isValid():
            continue
        for feature in triangle.getFeatures():
            fingerprint[tracker[feature]] += 1

    return fingerprint.tolist()
