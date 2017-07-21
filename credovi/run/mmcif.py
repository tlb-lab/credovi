"""
Prints all the current non-obsolete PDB codes found in the mmCIF database. The
PDB file mirror is not used directly because entries have to exist in the database
as well in order to identify entities, etc..
"""
from sqlalchemy import create_engine, text

from credovi import app
from credovi.schema import engine

def apply_offset(pdbs, offset):
    '''
    '''
    IN_OFFSET = False

    for pdb in pdbs:
        if pdb == offset: IN_OFFSET = True
        if IN_OFFSET: yield pdb

def apply_limit(pdbs, limit):
    '''
    '''
    # try to use limit as integer
    try:
        limit = int(limit)

        for i,pdb in enumerate(pdbs,1):
            if i <= limit: yield pdb
            else: break

    # limit is string
    except ValueError:
        for pdb in pdbs:
            if pdb == limit: break
            yield pdb

def get_current_pdbs(args):
    '''
    Every structure not in this list is most likely obsolete.
    '''
    statement = text("""
                       SELECT structure_id as pdb
                         FROM {mmcif}.entry e
                    LEFT JOIN pdb_dev.banned b ON e.structure_id = b.pdb
                        WHERE b.pdb IS NULL
                     ORDER BY 1
                     """.format(mmcif=app.config.get('database','mmcif')))

    engine.echo = False  # Forced echo off since SQLAlchemy prints to STDOUT and output gets mixed with the PDB list.
    result = engine.execute(statement).fetchall()

    pdbs = (row['pdb'] for row in result)

    if args.offset: pdbs = apply_offset(pdbs, args.offset.upper())
    if args.limit: pdbs = apply_limit(pdbs, args.limit)

    return pdbs

def do(controller):
    """
    """
    # get the controller command
    cmd = controller.command

    # get the command line arguments and options
    args = controller.pargs

    for pdb in get_current_pdbs(args):
        print pdb
