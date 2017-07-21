import os, sys
import csv
import gzip
import re
import logging

from os.path import getmtime
from datetime import datetime

from sqlalchemy import *
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY

# register new dialect to write tab-delimited files
csv.register_dialect('tabs', delimiter='\t')



try:
    from argparse import ArgumentParser
    from getpass import getpass
except ImportError:
    print >> sys.stderr, "This program requires Python version 2.7 or later"
    sys.exit(1)
else:
    parser = ArgumentParser()
    parser.add_argument('-u', '--user',   dest="dbuser", help="Username on database", default='bernardo')
    parser.add_argument('-H', '--host',   dest="dbhost", help="Location of the database server", default='bahamut')
    parser.add_argument('-p', '--port',   dest="dbport", type=int, help="Port of the database server", default=5432)
    parser.add_argument('-d', '--db',     dest="dbname", help="Name of the database", default='cryst')
    parser.add_argument('-s', '--schema', dest="schema", help="Schema to use on the database", default='pdb_dev')
    parser.add_argument('-r', '--resume', dest="skip_dump", action='store_true', help="Skip initial dump if file exists.", default=False)
    #parser.add_argument('-f', '--force',  dest="force", action='store_true', help="Force overwrites", default=False)
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Echo SQL commands", default=False)
    parser.add_argument('-y', '--yes', dest="confirm", action='store_true', help="Auto-confirm", default=False)
    parser.add_argument('-L', '--loglevel',dest="loglevel", action='store', choices=['debug','info','warn','error'],
                        help="Logging level to display (default: info)", default='info')

    opt = parser.parse_args()

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                    datefmt='%m/%d/%Y %I:%M:%S %p', level=getattr(logging, opt.loglevel.upper()))

PDB_MIRROR_DIR  = "/tlbnas/mirror/pdb/data/structures/divided/pdb"
BIOMT_DUMP      = '/tlbnas/temp/bahamut/biomt.%s' % opt.schema

engine      = create_engine(URL(drivername='postgresql+psycopg2', username=opt.dbuser,
                                password=getpass("DB Password:"), host=opt.dbhost,
                                port=opt.dbport, database=opt.dbname),
                            execution_options={'autocommit':True}, echo=opt.verbose)

connection      = engine.connect()
metadata        = MetaData(bind=engine)

def create_tables(metadata, schema):
    """
    """
    rawdata = Table('raw_biomt', metadata,
                    Column('pdb', String(4), nullable=False),
                    Column('assembly_serial', Integer),
                    Column('assembly_type', Text),
                    Column('determined_by', ARRAY(String)),
                    Column('pdb_chain_id', String(12)),
                    Column('operation_serial', Integer),
                    Column('rotation', ARRAY(Float)),
                    Column('translation', ARRAY(Float)),
                    schema=schema, prefixes=['unlogged'])

    Index('idx_rawdata_pdb', rawdata.c.pdb, rawdata.c.assembly_serial)

    biomt = Table('biomt', metadata,
                  Column('biomt_id', Integer, nullable=False),
                  Column('pdb', String(4), nullable=False),
                  Column('assembly_serial', Integer),
                  Column('assembly_size', Integer),
                  Column('assembly_type', Text),
                  Column('determined_by', ARRAY(String)),
                  Column('is_monomeric', Boolean(create_constraint=False), DefaultClause('false')),
                  Column('all_chains_at_identity', Boolean(create_constraint=False), DefaultClause('false')),
                  schema=schema)

    PrimaryKeyConstraint(biomt.c.biomt_id, deferrable=True, initially='deferred')
    Index('idx_biomt_assembly', biomt.c.pdb, biomt.c.assembly_serial, unique=True)

    biomt_ops = Table('biomt_ops', metadata,
                      Column('biomt_op_id', Integer, nullable=False),
                      Column('biomt_id', Integer, nullable=False),
                      Column('pdb_chain_id', String(12)),
                      Column('operation_serial', Integer),
                      Column('rotation', ARRAY(Float)),
                      Column('translation', ARRAY(Float)),
                      Column('is_at_identity', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      schema=schema)

    PrimaryKeyConstraint(biomt_ops.c.biomt_op_id, deferrable=True, initially='deferred')
    Index('idx_biomt_ops_biomt_id', biomt_ops.c.biomt_id, biomt_ops.c.pdb_chain_id)

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

def get_pdb_files():
    """
    """
    # get all PFB files from the mirror recursively
    for folder, subs, files in os.walk(PDB_MIRROR_DIR):
        for f in files:

            # make sure to only include structure files and nothing else
            if not f.endswith(".ent.gz"): continue

            yield f

def extract_remark350(path):
    """
    """
    # the lines of the REMARK 350 section will be appended to this list
    remark = []

    # the PDB structures are gzipped
    with gzip.open(path) as pdbfh:
        for line in pdbfh:
            if line.startswith("REMARK 350"):

                # remove newline and REMARK 350 at the beginning
                remark.append(line[11:].rstrip())

            # end parsing this file at the SEQRES section
            if line.startswith("SEQRES"):
                break

        # the first line will be an empty string!
        return remark[1:]

def extract_bio_oper(remark):
    """
    """
    # pattern th extract the chains that the transformation applies to
    pchn = re.compile("^(APPLY THE FOLLOWING TO|(?P<cntd>\s+AND)) CHAIN[S]?[:]\s?(?P<chains>.+)$")

    op = {}
    operation_id = 0

    for line in remark:

        # a new multimer starts here
        if line.startswith('BIOMOLECULE'):

            # there might be more than one space between the BIOMOLECULE: start
            # and the actual number
            try:
                biomol = int(line[12:])
            except ValueError:
                raise ValueError("cannot parse BIOMOLECULE: {}".format(line[12:]))

            op.update({biomol:{'chains':{}}})

            # attach by whom the quaternary structure was determined
            op[biomol]['determined_by'] = []

        # try to extract the determined biological unit
        elif line.startswith('AUTHOR DETERMINED BIOLOGICAL UNIT'):
            try:
                op[biomol]['determined_by'].append('author')
                op[biomol]['assembly'] = line[35:].lower()

                logging.debug("Author determined biological unit {}: {}"
                              .format(biomol, op[biomol]['assembly']))

            # the BIOMOLECULE: 1 line is missing from header
            except UnboundLocalError as e:
                raise RuntimeError("cannot parse REMARK 350! BIOMOLECULE serial "
                                   "is missing: {}".format(e.message))

        # try to extract the determined biological unit
        elif line.startswith('SOFTWARE DETERMINED QUATERNARY STRUCTURE'):
            try:
                op[biomol]['determined_by'].append('software')
                op[biomol]['assembly'] = line[42:].lower()

                logging.debug("Software determined biological unit {}: {}"
                              .format(biomol, op[biomol]['assembly']))

            # the BIOMOLECULE: 1 line is missing from header
            except UnboundLocalError as e:
                raise RuntimeError("cannot parse REMARK 350! BIOMOLECULE serial "
                                   "is missing: {}".format(e.message))

        # try to extract the list of chains the transformation applies to
        elif pchn.match(line):
            match = pchn.match(line)

            # extract the list of PDB chain identifiers
            _chains = [chain for chain in match.group('chains').split(',')]
            _chains = [chain.strip() for chain in _chains if chain.strip()]

            # add to the existing list of chains if this line continues with chains
            if match.group('cntd'): chains.extend(_chains)

            # otherwise create new list
            else: chains = _chains

            for chain in chains:
                if len(chain) == 1:
                    op[biomol]['chains'].update({chain:{}})
                else:
                    logging.warn("Thee chain {} is not valid and was not parsed "
                                 "properly!".format(chain))

            logging.debug("Operation for biomolecule {} found: chains {}"
                          .format(biomol, chains))

        # first row of rotations and translations
        elif line.startswith('  BIOMT1'):
            operation_id = int(line[10:12])

            # initialize matrices
            rot = map(float, line[13:42].split())
            trans = [float(line[46:])]

        elif line.startswith('  BIOMT2'):
            operation_id = int(line[10:12])
            rot.extend(map(float, line[13:42].split()))
            trans.append(float(line[46:]))

        elif line.startswith('  BIOMT3'):
            operation_id = int(line[10:12])
            rot.extend(map(float, line[13:42].split()))
            trans.append(float(line[46:]))

            # add matrices to dictionary for each chain and entity
            try:
                for chain in chains:
                    op[biomol]['chains'][chain].update({operation_id:(rot, trans)})
            except UnboundLocalError as e:
                raise RuntimeError("cannot assign operations to chains: {}"
                                   .format(e.message))

    return op

def main():
    """
    """

    if os.path.exists(BIOMT_DUMP) and opt.skip_dump:
        print >> sys.stderr, "Dump file %s exists. Reusing." % BIOMT_DUMP
        most_recent = datetime.fromtimestamp(getmtime(BIOMT_DUMP))
    else:
        fh = open(BIOMT_DUMP, 'w')
        writer = csv.writer(fh, dialect='tabs')

        logging.info("Commencing dump.")
        most_recent = datetime.fromtimestamp(0)

        for f in get_pdb_files():
            pdbid = f[3:7]
            path = os.path.join(PDB_MIRROR_DIR, pdbid[1:3], f)

            logging.info("starting to extract BIOMT of PDB entry %s.", pdbid.upper())

            if not os.path.isfile(path):
                logging.error("cannot extract BIOMT: path %s does not exist!", path)
                continue

            # extract the REMARK 350 lines from the PDB header
            remark = extract_remark350(path)

            # parse the section and turn it into a hierarchical form
            try:
                biomt = extract_bio_oper(remark)
            except StandardError, e:
                logging.error("Could not extract operations for PDB entry %s: %s", pdbid.upper(), e.message)
                continue

            # write transformations
            for assembly_serial, biomols in biomt.items():

                try: assembly = biomols.pop('assembly')
                except KeyError: assembly = '\N'

                try:
                    determined_by = biomols.pop('determined_by')
                    determined_by = '{' + ','.join(determined_by) + '}'
                except KeyError:
                    determined_by = '\N'

                for pdb_chain_id, operation in biomols['chains'].items():
                    for operation_serial, (rx,tx) in operation.items():
                        rx = '{' + ','.join(map(str,rx)) + '}'
                        tx = '{' + ','.join(map(str,tx)) + '}'

                        fields = [pdbid.upper(), assembly_serial, assembly, determined_by,
                                  pdb_chain_id, operation_serial, rx, tx]

                        writer.writerow(fields)

            pdb_mtime = datetime.fromtimestamp(getmtime(path))
            if pdb_mtime > most_recent:
                most_recent = pdb_mtime
            fh.flush()
        fh.close()

    if opt.confirm or raw_input("Tables in schema {} are about be dropped and recreated. " \
                                "Are you sure you want to continue? [y/N] ".format(opt.schema)).strip().lower() == 'y':
        logging.info("Creating tables...")
        create_tables(metadata, opt.schema)
    else:
        print >> sys.stderr, "Aborting."
        sys.exit(1)

    logging.info("Copying dump file to table...")

    connection.execute("COPY {}.raw_biomt FROM '{}'".format(opt.schema, BIOMT_DUMP))

    logging.info("Inserting into {}.biomt...".format(opt.schema))

    connection.execute("""INSERT INTO {schema}.biomt(pdb, assembly_serial, assembly_type, determined_by)
                            SELECT DISTINCT pdb, assembly_serial, assembly_type, determined_by
                              FROM {schema}.raw_biomt
                          ORDER BY 1,2;""".format(schema=opt.schema))

    logging.info("Inserting into %s.biomt_ops...", opt.schema)

    connection.execute("""INSERT INTO {schema}.biomt_ops(biomt_id, pdb_chain_id, operation_serial, rotation, translation)
                            SELECT biomt_id, pdb_chain_id, operation_serial, rotation, translation
                              FROM {schema}.biomt
                              JOIN {schema}.raw_biomt USING(pdb, assembly_serial)
                          ORDER BY 1,2,3;""".format(schema=opt.schema))

    connection.execute("COMMENT ON TABLE {schema}.biomt IS 'Most recent data: {timestamp}'"\
                       .format(schema=opt.schema, timestamp=str(most_recent)) )
    connection.execute("COMMENT ON TABLE {schema}.biomt_ops IS 'Most recent data: {timestamp}'"\
                       .format(schema=opt.schema, timestamp=str(most_recent)) )


    connection.execute("""UPDATE {schema}.biomt_ops
                          SET is_at_identity = true
                          WHERE rotation[1] = 1 and rotation[5] = 1 and rotation[9] = 1 and sum(rotation) = 3
                              and translation[1] = 0 and translation[2] = 0 and translation[3] = 0;""".format(schema=opt.schema))

    connection.execute(""" UPDATE {schema}.biomt b
                            SET all_chains_at_identity = sq.all_chains_at_identity
                           FROM (
                                   SELECT biomt_id, bool_and(is_at_identity) as all_chains_at_identity
                                     FROM {schema}.biomt_ops
                                 GROUP BY biomt_id
                                ) sq
                          WHERE sq.biomt_id = b.biomt_id;""".format(schema=opt.schema))

    connection.execute(""" UPDATE {schema}.biomt b
                            SET assembly_size = sq.assembly_size
                           FROM (
                                   SELECT biomt_id, count(operation_serial) as assembly_size
                                     FROM {schema}.biomt_ops
                                 GROUP BY biomt_id
                                ) sq
                          WHERE sq.biomt_id = b.biomt_id""".format(schema=opt.schema))

    connection.execute(""" UPDATE {schema}.biomt b
                            SET is_monomeric = true
                           FROM (
                                   SELECT pdb
                                     FROM {schema}.biomt
                                 GROUP BY pdb
                                   HAVING bool_and(all_chains_at_identity) = true
                                          AND max(assembly_size) = 1
                                ) sq
                          WHERE sq.pdb = b.pdb""".format(schema=opt.schema))

main()
