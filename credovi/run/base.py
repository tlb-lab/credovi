"""
This module is used to query the current configuration settings.
"""

import os
from collections import Iterable

from credovi import app
from credovi.util import filesystem as fs

def printvalue(args, value):
    """
    """
    if args.sort: value.sort()
        
    # if the item is iterable, print each item on a new line
    if args.single and isinstance(value, Iterable):
        for item in value:
            print item
    
    # otherwise print the complete object
    else: print value

def conf(args):
    """
    Print sections and keys of the current CREDO configuration.
    """
    section, keys = args.section, args.keys
        
    if section:
        if app.config.has_section(section):            

            # print the specified keys of the section
            if keys:
                for key in keys.split(','):
                    if app.config.has_key(section, key):
                        value = app.config.get(section, key)
                    
                        printvalue(args, value)
                    
                    # key does not exist in the specified section
                    else:
                        app.log.warn("key [{0}] does not exist in section "
                                     "[{1}]".format(key, section))
        
            # print all keys in the section
            else:
                for key in app.config.keys(section):
                    value = app.config.get(section, key)

                    printvalue(args, value)
        
        # the specified section does not exist
        else:
            app.log.error("section [{0}] does not exist in the current "
                          "configuration.".format(section))

def data(args):
    """
    Work with the local CREDO data directory.
    """
    # fetch all PDB entries from the CREDO structures table and compare with
    # assembly sets in CREDO data directory
    if args.list_not_in_db:

        # finally print PDB codes: one per line, parallel friendly
        for pdb in fs.get_pdbs_not_in_credodb(): print pdb

def do(controller):
    """
    """
    # get the controller command
    cmd = controller.command
    
    # get the command line arguments and options
    args = controller.pargs  
    
    # print sections and keys of the current CREDO configuration
    if cmd == 'conf': conf(args)
    
    # work with the CREDO data directory
    elif cmd == 'data': data(args)