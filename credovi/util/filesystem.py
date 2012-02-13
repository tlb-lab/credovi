"""
This module contains functions to interact with the filesystem, e.g. the CREDO
data directory.
"""

import os

from credovi import app
from credovi.structbio import db

def get_credo_data_dir_pdbs():
    """
    Returns a frozenset containing all the PDB entries that can be found in the
    CREDO data directory. 
    """
    credo_data_dir = app.config['directories']['credo_data']
    
    pdbs = set()
    
    # CREDO data directory is divided so we need to use os.walk to get into the
    # subdirectories as well
    for pdbgroupdir in os.listdir(credo_data_dir):      
        
        # directory is not two characters long
        if len(pdbgroupdir) != 2:
            app.log.warn("PDB group directory {} is not valid!".format(pdbgroupdir))
            continue
        
        # we need the full path again for os.listdir()
        pdbgroupdir = os.path.join(credo_data_dir, pdbgroupdir)            
        
        # each directory here is the lower case PDB code
        for pdbdir in os.listdir(pdbgroupdir): pdbs.add(pdbdir.upper())

    # debug how many PDB entry we found
    app.log.debug("{} PDB entries were found in CREDO data directory."
                  .format(len(pdbs)))    
    
    return pdbs

def get_pdbs_not_in_credodb():
    """
    Returns a list of all the PDB entries in the CREDO data directory that are
    not part of the CREDO database yet.
    """
    # get the PDB entries in the CREDO data directory as set
    pdbs = get_credo_data_dir_pdbs()
    
    # iterate through all the PDB entries in the CREDO database and remove
    # them from the set
    for pdb in db.get_credo_pdbs():
        if pdb in pdbs:
            pdbs.discard(pdb)
        
        # PDB entry is in the database but not in the CREDO data directory
        else:
            app.log.warn("PDB entry {} is in the database but not in the CREDO "
                         "data directory!".format(pdb))

    # info how many PDB entries are left after intersecting the ones from the db
    app.log.debug("{} PDB entries in the CREDO data directory are not in "
                  "the database yet.".format(len(pdbs)))
    
    return sorted(pdbs)