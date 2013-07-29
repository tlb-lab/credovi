import os
import csv
import gzip
import re
import logging

from sqlalchemy import *
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY

# register new dialect to write tab-delimited files
csv.register_dialect('tabs', delimiter='\t')

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                    datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

PDB_MIRROR_DIR  = "/tlbnas/mirror/pdb/data/structures/divided/pdb"
BIOMT_DUMP      = '/tlbnas/temp/bahamut/biomt.pdb'

engine          = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                    password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                    port=5432, database='cryst'),
                                execution_options={'autocommit':True}, echo=False)
connection      = engine.connect()
metadata        = MetaData(bind=engine)

def create_tables():
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
                    schema='pdb', prefixes=['unlogged'])

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
                  schema='pdb')

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
                      schema='pdb')

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
    fh = open(BIOMT_DUMP, 'w')
    writer = csv.writer(fh, dialect='tabs')

    for f in get_pdb_files():
        pdb = f[3:7]
        path = os.path.join(PDB_MIRROR_DIR, pdb[1:3], f)

        logging.info("starting to extract BIOMT of PDB entry {}."
                     .format(pdb.upper()))

        if not os.path.isfile(path):
            logging.error("cannot extract BIOMT: path {} does not exist!"
                          .format(path))
            continue

        # extract the REMARK 350 lines from the PDB header
        remark = extract_remark350(path)

        # parse the section and turn it into a hierarchical form
        try:
            biomt = extract_bio_oper(remark)
        except (ValueError, RuntimeError) as e:
            logging.error("cannot extract operations for PDB entry {}: {}"
                          .format(pdb.upper(), e.message))
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

                    fields = [pdb.upper(), assembly_serial, assembly, determined_by,
                              pdb_chain_id, operation_serial, rx, tx]

                    writer.writerow(fields)

        fh.flush()
    fh.close()

    logging.info("Creating tables...")

    create_tables()

    logging.info("Copying dump file to table...")

    connection.execute("COPY pdb.raw_biomt FROM '{}'".format(BIOMT_DUMP))

    logging.info("Inserting into pdb.biomt...")

    connection.execute("""
                            INSERT INTO pdb.biomt(pdb, assembly_serial, assembly_type, determined_by)
                            SELECT DISTINCT pdb, assembly_serial, assembly_type, determined_by
                              FROM pdb.raw_biomt
                          ORDER BY 1,2;""")

    logging.info("Inserting into pdb.biomt_ops...")

    connection.execute("""
                            INSERT INTO pdb.biomt_ops(biomt_id, pdb_chain_id, operation_serial, rotation, translation)
                            SELECT biomt_id, pdb_chain_id, operation_serial, rotation, translation
                              FROM pdb.biomt
                              JOIN pdb.raw_biomt USING(pdb, assembly_serial)
                          ORDER BY 1,2,3;
                       """)

    connection.execute("""
                       UPDATE pdb.biomt_ops
                          SET is_at_identity = true
                        WHERE rotation[1] = 1 and rotation[5] = 1 and rotation[9] = 1 and sum(rotation) = 3
                              and translation[1] = 0 and translation[2] = 0 and translation[3] = 0;
                       """)

    connection.execute("""
                         UPDATE pdb.biomt b
                            SET all_chains_at_identity = sq.all_chains_at_identity
                           FROM (
                                   SELECT biomt_id, bool_and(is_at_identity) as all_chains_at_identity
                                     FROM pdb.biomt_ops
                                 GROUP BY biomt_id
                                ) sq
                          WHERE sq.biomt_id = b.biomt_id;
                       """)

    connection.execute("""
                         UPDATE pdb.biomt b
                            SET assembly_size = sq.assembly_size
                           FROM (
                                   SELECT biomt_id, count(operation_serial) as assembly_size
                                     FROM pdb.biomt_ops
                                 GROUP BY biomt_id
                                ) sq
                          WHERE sq.biomt_id = b.biomt_id
                       """)

    connection.execute("""
                         UPDATE pdb.biomt b
                            SET is_monomeric = true
                           FROM (
                                   SELECT pdb
                                     FROM pdb.biomt
                                 GROUP BY pdb
                                   HAVING bool_and(all_chains_at_identity) = true
                                          AND max(assembly_size) = 1
                                ) sq
                          WHERE sq.pdb = b.pdb
                       """)

main()
