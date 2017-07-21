"""
"""
import os
import string
from glob import iglob
from itertools import groupby

from sqlalchemy.schema import CreateIndex, CreateTable

from credovi import app
from credovi.schema import metadata, schema
from credovi.util import filesystem as fs
from credovi.util import postgresql as pg

def get_credo_files(pdbs):
    """
    """
    # make sure that all PDB codes are lower case
    pdbs = map(string.lower, pdbs)

    for pdb in pdbs:
        pdbdir = os.path.join(app.config['directories']['credo_data'], pdb[1:3], pdb)

        # get the .credo files in the CREDO data directory for this PDB entry
        for path in iglob(os.path.join(pdbdir, '*.credo')):
            yield path

def copy(args):
    """
    Copy the data files of all PDB entries that are not in CREDO yet.
    """
    # copy new CREDO data into the raw tables of the database
    if args.new:

        # get all the PDB codes of entries not in the CREDO database
        pdbs = fs.get_pdbs_not_in_credodb()

    else: pdbs = []

    # get an iterator over the .credo files of all PDB entries
    credofiles = get_credo_files(pdbs)

    # key function to sort & group credo files by name
    key = lambda path: os.path.split(path)[1]

    # generator cannot be sorted / should fine a different solution here
    credofiles = sorted(credofiles, key=key)

    for filename, groupiter in groupby(credofiles, key=key):

        # choose the tablename depending on the file type
        if filename == 'ligands.credo': table = 'credo.raw_ligands'
        elif filename == 'atoms.credo': table = 'credo.raw_atoms'
        elif filename == 'contacts.credo': table = 'credo.raw_contacts'
        elif filename == 'aromaticrings.credo': table = 'credo.raw_aromaticrings'
        else:
            app.log.warn("cannot copy uknown file type {} into database!"
                         .format(filename))
            continue

        # do not use multiprocessing with VirtualBox
        pg.copy(groupiter, table=table, processes=1)

def create(args, tablenames=None):
    """
    """
    if not tablenames:

        if args.raw or args.core:
            for table in metadata.tables.values():

                # only create the 'raw' tables used for loading data
                if args.raw and table.name.startswith('raw'):
                    table.drop(checkfirst=args.checkfirst)
                    table.create(checkfirst=args.checkfirst)

                # core tables
                elif args.core and table.name.startswith('raw'):
                    metadata.remove(table)

            if args.core:
                metadata.drop_all(checkfirst=args.checkfirst)
                metadata.create_all(checkfirst=args.checkfirst)

        # create all tables
        else:
            if args.sure:
                app.log.info("creating all elements defined in the current schema!")

                metadata.drop_all(checkfirst=args.checkfirst)
                metadata.create_all(checkfirst=args.checkfirst)
            else:
                app.log.error("you must specify the --sure option before attempting "
                              "to drop and create the whole schema.")
                app.close()

    else:
        for tablename in tablenames:
            if tablename in metadata.tables:
                table = metadata.tables[tablename]
                table.drop(checkfirst=args.checkfirst)
                table.create(checkfirst=args.checkfirst)

            else:
                #current_names = metadata.tables.keys()
                app.log.warn("cannot create table {0}: not defined in the current schema.".format(tablename))

def drop(args, tablenames=None):
    """
    """
    if not tablenames:
        app.log.debug("dropping all elements defined in the current schema!")

        metadata.drop_all(checkfirst=args.checkfirst)

    else:
        for tablename in tablenames:
            if tablename in metadata.tables:
                table = metadata.tables[tablename]
                table.drop(checkfirst=args.checkfirst)

            else:
                app.log.error("cannot drop table {0}: not defined in the current schema.")

def dump(args, tablenames=None):
    """
    """
    # complete schema
    if not tablenames:
        app.log.debug("dumping metadata of all tables in the CREDO schema.")

        for table in metadata.tables.values():
            print CreateTable(table, on='postgresql', bind=metadata.bind)

            for index in table.indexes:
                print CreateIndex(index, on='postgresql', bind=metadata.bind)

    # only specific tables
    else:

        # check if the given tables are are defined in CREDO
        for tablename in tablenames:
            if tablename not in metadata.tables:
                app.log.fatal("table {0} is not defined in CREDO".format(tablename))
                app.close()

            else:
                table = metadata.tables[tablename]
                print CreateTable(table, on='postgresql', bind=metadata.bind)

                for index in table.indexes:
                    print CreateIndex(index, on='postgresql', bind=metadata.bind)

def do(controller):
    """
    Here is the application logic for dealing with the CREDO database schema. The
    command line interface (CLI) controller is parsed from the main script.
    """
    # get the controller command
    cmd = controller.command

    # get the command line arguments and options
    args = controller.pargs

    # override the SQL echo setting if specified on command line
    if args.echo: metadata.bind.echo = args.echo

    # get all tablenames that were given on the command line as list
    if args.tables:
        tablenames = ['.'.join((schema, name)) if schema and not name.startswith(schema) else name
                       for name in args.tables.split(',')]
        if args.add_depend:
            tablenames.extend([tab for tab in metadata.tables
                               if tab not in tablenames and any(tab.startswith(n) for n in tablenames)])
    else:
        tablenames = None

    # no command was given, exit
    if not cmd: pass

    #
    elif cmd == 'copy': copy(args)

    # create elements in the CREDO database schema
    elif cmd == 'create': create(args, tablenames)

    # drop the database schema
    elif cmd == 'drop': drop(args, tablenames)

    # dump the currently defined CREDO database schema
    elif cmd == 'dump': dump(args, tablenames)


    # truncate tables in the database
    elif cmd == 'truncate': pass
