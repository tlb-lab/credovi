import sys, os
import gzip
import json
import csv
import logging
import multiprocessing as mp
from getpass import getpass
from copy import copy
from os.path import getmtime, exists, getsize, isfile
from datetime import datetime

###############################################################################################
####  Along with create-biomt-from-header.py, this script populates the 'pdb[_dev]' schema for CREDO.
###############################################################################################
try:
    from lxml import etree # Changed xml lib for lxml
except ImportError: # Python 2.5, C version
    print >> sys.stderr, "lxml could not be imported. Automatic namespace definition not possible."
    import xml.etree.cElementTree as etree

from sqlalchemy import *
from sqlalchemy.event import listen
#from sqlalchemy.schema import CreateIndex
from sqlalchemy.exc import ProgrammingError, DatabaseError, DataError, SQLAlchemyError
from sqlalchemy.engine.url import URL
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

# CONFIG      = json.loads(open('credovi/config/project.json').read()) # unused, non-existent, config.json?
CREDO_CONF  = json.loads(open('credovi/config/credo.json').read())
SIFTS_DIR   = '/tlbnas/mirror/pdbe/sifts/xml'
#NS         = '{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}'               # Deprecated because not all SIFT files are consistent.
#NS         = '{http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd}'   # Now using dynamic automatic definition (requires lxml library)

if __name__ == '__main__':

    try:
        from argparse import ArgumentParser
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
        parser.add_argument('-m', '--mmcif',  dest="mmcif_schema", help="Schema on the database where MMCIF data is stored", default='mmcif_dev')
        parser.add_argument('-r', '--resume', dest="skip_dump", action='store_true', help="Skip initial dump if file exists.", default=False)
        parser.add_argument('--skip-load', dest="skip_load", action='store_true', help="Skip loading of dump file (i.e. only create indices + frag/lig tables).", default=False)
        #parser.add_argument('-f', '--force',  dest="force", action='store_true', help="Force overwrites", default=False)
        parser.add_argument('-c', '--cit-only', dest="citations_only", action='store_true', help="Only process citations", default=False)
        parser.add_argument('-T', '--test', dest="test", nargs='?', action='store', const=True, help="Process test set only", default=False)
        parser.add_argument('-y', '--yes', dest="confirm", action='store_true', help="Auto-confirm", default=False)
        parser.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Echo SQL commands", default=False)
        parser.add_argument('-L', '--loglevel',dest="loglevel", action='store', choices=['debug','info','warn','error'],
                            help="Logging level to display (default: info)", default='info')
    
        opt = parser.parse_args()
    
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=getattr(logging, opt.loglevel.upper()))
    
    engine      = create_engine(URL(drivername='postgresql+psycopg2', username=opt.dbuser,
                                    password=getpass("DB Password for %s:" % opt.dbuser), host=opt.dbhost,
                                    port=opt.dbport, database=opt.dbname),
                                execution_options={'autocommit':True}, echo=opt.verbose)
    connection  = engine.connect()
    metadata    = MetaData(bind=connection)
    
    
    try:
        if opt.test == True:
            PDB_TEST_SET = frozenset(CREDO_CONF['test sets']['large'])
        elif opt.test:
            PDB_TEST_SET = opt.test.split(',')
        else:
            PDB_TEST_SET = [] #['PDB TEST SET'])
    except KeyError, err:
        logging.error("Test set key not found on %s: %s.", CREDO_CONF, err)
        sys.exit(1)
    
    SCHEMA = opt.schema
    MMCIF  = opt.mmcif_schema
else:
    SCHEMA = 'pdb_dev'
    MMCIF  = 'mmcif_dev'
    
    logging.basicConfig(format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=getattr(logging, os.environ.get('loglevel', 'info').upper()))
    
    engine      = create_engine(URL(drivername='postgresql+psycopg2', username=raw_input('DB Username: ').strip(),
                                    password=getpass("DB Password:"), host='bahamut',
                                    port=5432, database='cryst'),
                                execution_options={'autocommit':True}, echo=False)
    connection  = engine.connect()
    metadata    = MetaData(bind=connection)


RES_MAP     = '/tlbnas/temp/bahamut/res_map.%s' % SCHEMA
MAP_REGIONS = '/tlbnas/temp/bahamut/map_regions.%s' % SCHEMA
CITATIONS   = '/tlbnas/temp/bahamut/citations.%s' % SCHEMA

def reflect_tables(metadata):
    res_map     = Table('res_map',     metadata, schema=SCHEMA, autoload=True)
    map_regions = Table('map_regions', metadata, schema=SCHEMA, autoload=True)
    citations   = Table('citations',   metadata, schema=SCHEMA, autoload=True)
    
    logging.info("Tables res_map, map_regions and citations successfully reflected.")
    return res_map, map_regions, citations

def create_tables(metadata=metadata):
    '''
    '''
    res_map = Table('res_map', metadata,
                    Column('res_map_id', Integer, primary_key=True),
                    Column('entry', String(4), nullable=False),
                    Column('entity_id', String(3), nullable=False),
                    Column('sifts_res_num', Integer, nullable=False),
                    Column('sifts_res_name', String(3)),
                    Column('pdb', String(4)),
                    Column('pdb_chain_id', String(4)),
                    Column('pdb_res_num', Integer),
                    Column('pdb_ins_code', String(1)),
                    Column('pdb_res_name', String(3)),
                    Column('uniprot', String(10)),
                    Column('uniprot_res_num', Integer),
                    Column('uniprot_res_name', String(1)),
                    Column('sstruct', String(1)),
                    Column('sstruct_serial', Integer),
                    Column('is_observed', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    Column('is_modified', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    Column('is_conflicting', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                    schema=SCHEMA)

    map_regions = Table('map_regions', metadata,
                        Column('map_region_id', Integer, primary_key=True),
                        Column('pdb', String(4), nullable=False),
                        Column('pdb_chain_id', String(4), nullable=False),
                        Column('sifts_res_num_start', Integer, nullable=False),
                        Column('sifts_res_num_end', Integer, nullable=False),
                        Column('db_source', String(12)),
                        Column('db_accession_id', String(20)),
                        Column('db_version', String(12)),
                        Column('db_coord_sys', String(12)),
                        schema=SCHEMA)

    disordered_regions = Table('disordered_regions', metadata,
                               Column('map_region_id', Integer, primary_key=True),
                               Column('pdb', String(4), nullable=False),
                               Column('pdb_chain_id', String(4), nullable=False),
                               Column('disordered_region_serial', Integer, nullable=False),
                               Column('region_seq', Text),
                               Column('region_length', Integer, nullable=False),
                               Column('region_start', Integer, nullable=False),
                               Column('region_end', Integer, nullable=False),
                               Column('is_nterminus', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                               Column('is_cterminus', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                               schema=SCHEMA)

    Index('idx_disordered_regions_pdb', disordered_regions.c.pdb, disordered_regions.c.pdb_chain_id)

    pdb_prot_fragments = Table('pdb_prot_fragments', metadata,
                               Column('pdb_prot_fragment_id', Integer, primary_key=True),
                               Column('pdb', String(4), nullable=False),
                               Column('pdb_chain_id', String(4), nullable=False),
                               Column('sstruct_serial', Integer),
                               Column('sstruct', String(1)),
                               Column('fragment_size', Integer),
                               Column('fragment_seq', Text),
                               Column('fragment_nterm', Integer),
                               Column('fragment_cterm', Integer),
                               schema=SCHEMA)

    Index('idx_pdb_prot_fragments_pdb', pdb_prot_fragments.c.pdb, pdb_prot_fragments.c.pdb_chain_id, pdb_prot_fragments.c.sstruct_serial, unique=True)

    pdb_prot_fragment_to_residue = Table('pdb_prot_fragment_to_residue', metadata,
                                         Column('pdb_prot_fragment_id', Integer, primary_key=True, autoincrement=False),
                                         Column('res_map_id', Integer, primary_key=True, autoincrement=False),
                                         schema=SCHEMA)

    Index('idx_pdb_prot_fragment_to_residue_res_map_id', pdb_prot_fragment_to_residue.c.res_map_id, unique=True)

    ## uniprot_residues don't seem to be used anywhere on the current code, hence disabled.
    # uniprot_residues = Table('uniprot_residues', metadata,
                             # Column('uniprot_residue_id', Integer, primary_key=True),
                             # Column('uniprot', String(6)),
                             # Column('uniprot_res_num', Integer),
                             # Column('uniprot_res_name', String(1)),
                             # schema=SCHEMA)


    peptide_ligands = Table('peptide_ligands', metadata,
                             Column('pdb', String(4), primary_key=True),
                             Column('pdb_chain_id', String(4), primary_key=True),
                             Column('description', Text),
                             Column('fragment', Text),
                             Column('seq', Text),
                             schema=SCHEMA)

    ligands = Table('ligands', metadata,
                    Column('pdb', String(4), primary_key=True),
                    Column('pdb_chain_id', String(4), primary_key=True),
                    Column('res_num', Integer, primary_key=True),
                    Column('ins_code', String(1), default=' ', primary_key=True),
                    Column('het_id', String(3), nullable=False, primary_key=True),
                    Column('name', Text),
                    Column('is_observed', Boolean, default=True),  # seems pointless (always true?)
                    Column('is_valid', Boolean,    default=None),  # seems pointless (always null?)
                    schema=SCHEMA)
    Index('idx_ligands_het_id', ligands.c.het_id) #, unique=True)

    citations = Table('citations', metadata,
                      Column('pubmed_id', Integer, primary_key=True),
                      Column('source', Text),
                      Column('title', Text),
                      Column('abstract', Text),
                      schema=SCHEMA)


    # MAKE THESE TABLES PUBLICLY ACCESSIBLE
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.res_map TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.map_regions TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.disordered_regions TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.pdb_prot_fragments TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.pdb_prot_fragment_to_residue TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.peptide_ligands TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.ligands TO public".format(SCHEMA)))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE {}.citations TO public".format(SCHEMA)))

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

    return res_map, map_regions, citations #, peptide_ligands, ligands

def update_sstruct_serials(connection=connection):
    '''
    '''
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM {schema}.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        UPDATE  {schema}.res_map m
                        SET     sstruct_serial = sq.sstruct_serial
                        FROM    (
                                WITH T1 AS
                                (
                                    SELECT  res_map_id, entry, entity_id, sifts_res_num, sstruct,
                                            LAG(sstruct,1,sstruct)
                                                OVER (ORDER BY entry, entity_id, sifts_res_num)
                                                IS DISTINCT FROM sstruct AS changes
                                    FROM {schema}.res_map WHERE pdb = $1
                                )
                                SELECT  res_map_id, entry, entity_id, sifts_res_num, sstruct,
                                        SUM(changes::int) OVER (PARTITION BY entry, entity_id
                                                                ORDER BY entry, entity_id, sifts_res_num) + 1 as sstruct_serial
                                FROM    T1
                                ) sq
                        WHERE   sq.res_map_id = m.res_map_id
                        ' USING pdbcode;
                    END LOOP;
                END$$;
                """.format(schema=SCHEMA)

    logging.info("Updating sstruct serials...")
    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

def insert_disordered_regions(connection=connection, timestamp=None):
    '''
    '''
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM {schema}.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        INSERT INTO {schema}.disordered_regions(pdb, pdb_chain_id, disordered_region_serial, region_seq, region_length, region_start, region_end, is_nterminus, is_cterminus)
                        WITH        W1 AS
                                    (
                                    SELECT  pdb, pdb_chain_id, pdb_res_num, bio.to_one_letter_code(pdb_res_name) as pdb_res_name, is_observed,
                                            LAG(is_observed, 1, true)
                                                OVER (PARTITION BY pdb, pdb_chain_id
                                                    ORDER BY pdb, pdb_chain_id, pdb_res_num)
                                                IS DISTINCT FROM is_observed AS changes,
                                            CASE
                                                WHEN LAG(pdb_res_num) OVER (PARTITION BY pdb, pdb_chain_id ORDER BY pdb, pdb_chain_id, pdb_res_num) IS NULL THEN 1
                                                ELSE 0
                                            END AS is_nterminus,
                                            CASE
                                                WHEN LEAD(pdb_res_num) OVER (PARTITION BY pdb, pdb_chain_id ORDER BY pdb, pdb_chain_id, pdb_res_num) IS NULL THEN 1
                                                ELSE 0
                                            END AS is_cterminus
                                    FROM    {schema}.res_map
                                    WHERE   pdb = $1
                                    ),
                                    W2 AS
                                    (
                                    SELECT  pdb, pdb_chain_id, pdb_res_num, pdb_res_name, is_observed,
                                            SUM(changes::int) OVER (PARTITION BY pdb, pdb_chain_id
                                                                    ORDER BY pdb, pdb_chain_id, pdb_res_num) as disordered_region_serial,
                                            is_nterminus, is_cterminus
                                    FROM    W1
                                    WHERE   is_observed = false
                                    )
                        SELECT      pdb, pdb_chain_id, disordered_region_serial,
                                    array_to_string(array_agg(pdb_res_name ORDER BY pdb_res_num),'''') as region_seq,
                                    COUNT(pdb_res_num) as region_length,
                                    min(pdb_res_num) as region_start, max(pdb_res_num) as region_end,
                                    max(is_nterminus)::boolean AS is_nterminus, max(is_cterminus)::boolean AS is_cterminus
                        FROM        W2
                        WHERE       is_observed = false
                        GROUP BY    pdb, pdb_chain_id, disordered_region_serial
                        ORDER BY    1,2,3
                        ' USING pdbcode;

                    END LOOP;
                END$$;
                """.format(schema=SCHEMA)

    logging.info("Inserting disordered regions...")
    transaction = connection.begin()
    connection.execute(statement)
    if timestamp:
        connection.execute("COMMENT ON TABLE {schema}.disordered_regions IS 'Most recent data: {timestamp}'".format(
            schema=SCHEMA, timestamp=str(timestamp)) )
    transaction.commit()

def insert_prot_fragments(connection=connection, timestamp=None):
    '''
    '''
    # CREATE UNIQUE REPRESENTATION FOR PROTEIN FRAGMENTS IN THE PDB
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM {schema}.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        INSERT      INTO {schema}.pdb_prot_fragments(pdb, pdb_chain_id, sstruct_serial, sstruct, fragment_size, fragment_seq, fragment_nterm, fragment_cterm)
                        SELECT      pdb, pdb_chain_id, sstruct_serial, sstruct,
                                    LENGTH(ARRAY_TO_STRING(ARRAY_AGG(bio.to_one_letter_code(sifts_res_name)),'''')) AS fragment_size,
                                    ARRAY_TO_STRING(array_agg(bio.to_one_letter_code(sifts_res_name) ORDER BY pdb_res_num),'''') AS seq,
                                    LAG(m.sstruct_serial) OVER(PARTITION BY pdb) AS nterminal,
                                    LEAD(m.sstruct_serial) OVER(PARTITION BY pdb) AS cterminal
                        FROM        {schema}.res_map m
                        WHERE       m.pdb = $1
                        GROUP BY    pdb, pdb_chain_id, sstruct_serial, sstruct
                        ORDER BY    pdb, pdb_chain_id, sstruct_serial
                        ' USING pdbcode;

                    END LOOP;
                END$$;
                """.format(schema=SCHEMA)

    logging.info("Inserting protein fragments...")
    transaction = connection.begin()
    connection.execute(statement)
    if timestamp:
        connection.execute("COMMENT ON TABLE {schema}.pdb_prot_fragments IS 'Most recent data: {timestamp}'".format(
            schema=SCHEMA, timestamp=str(timestamp)) )
    transaction.commit()

    # CREATE A MAPPING BETWEEN PDB PROTEIN FRAGMENTS AND THE RESIDUES FROM RESMAP
    statement = """
                INSERT      INTO {schema}.pdb_prot_fragment_to_residue
                SELECT      f.pdb_prot_fragment_id, m.res_map_id
                FROM        {schema}.pdb_prot_fragments f
                JOIN        {schema}.res_map m
                            ON m.pdb = f.pdb
                            AND m.pdb_chain_id = f.pdb_chain_id
                            AND m.sstruct_serial = f.sstruct_serial
                ORDER BY    f.pdb_prot_fragment_id, m.res_map_id
                """.format(schema=SCHEMA)

    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

def dump_citations(in_q, exc_l, connection=connection):
    import urllib, urllib2
    import Queue
    from urllib2 import HTTPError, URLError
    from sqlalchemy.sql import select, and_
    from time import sleep

    citation_fh = open(CITATIONS, 'w')
    csv.register_dialect('tabs', delimiter='\t')
    citationwriter = csv.DictWriter(citation_fh, restval='\N', extrasaction='ignore',
                                    dialect='tabs',  fieldnames = ['pubmed_id', 'source', 'title', 'abstract'])

    citations = Table('citation', metadata, schema=MMCIF, autoload=True)

    seen_pmid = set()
    pdbs_completed, pdbs_failed, pdbs_duped = 0, 0, 0
    conn_attempts = 0
    while True:
        try:
            pdbid = in_q.get(timeout=1800)
            #mp_logger.debug("Got arguments: %s", args)
            if pdbid is None:
                #mp_logger.debug("Found sentinel. Breaking.")
                break

            s = select([citations]).where(and_(citations.c.structure_id == pdbid, citations.c.id == 'primary'))
            row = connection.execute(s).fetchone()

            if not row:
                logging.debug('Could not find citation for PDB %s on MMCIF database', pdbid)
                pdbs_failed += 1
                exc_l.append((pdbid, "No citation"))
                continue

            rowdict = {}
            rowdict['pubmed_id'] = pmid = row['pdbx_database_id_pubmed']
            rowdict['title']     = row['title']
            if not pmid or int(pmid) <= 0:
                logging.debug('PubMed ID for PDB %s is not available on MMCIF database (%s)', pdbid, pmid)
                pdbs_failed += 1
                exc_l.append((pdbid, "No PubMed ID"))
                continue
            elif int(pmid) in seen_pmid:
                logging.debug('PubMed ID %s for PDB %s already processed.', pmid, pdbid)
                pdbs_duped += 1
                continue

            pmc_info_url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/resulttype=core&query=ext_id:%s&format=json" % (pmid)
            conn_pending = True
            while conn_pending or conn_attempts <= 3:
                try:
                    pmc_info_f = urllib2.urlopen(pmc_info_url, timeout=60)
                except HTTPError, err:
                    logging.warning("Could not retrieve PMC info for PDB %s, PMID %s", pdbid, pmid, exc_info=True)
                    pdbs_failed += 1
                    exc_l.append((pdbid, err.reason))
                    conn_pending = False
                except URLError, err:
                    logging.warning("Error retrieving PMC info for PDB %s, PMID %s: %s. Cooling off...",
                                    pdbid, pmid, err.reason)
                    sleep(30)
                    conn_attempts += 1
                else:
                    pmc_json = json.load(pmc_info_f)
                    break
            else:
                if conn_attempts > 3:
                    logging.error("Failed to connect to EBI server %d times. Aborting.", conn_attempts)
                    break
                conn_attempts = 0
                continue

            try:
                result   = pmc_json['resultList']['result'][0]
                if 'issue' in result['journalInfo']:
                    issue = '(%s)' % result['journalInfo']['issue']
                else:
                    issue = ''
                date     = result['journalInfo']['dateOfPublication']
                if 'abstractText' in result:
                    rowdict['abstract'] = result['abstractText']
                else:
                    logging.info("Missing abstract text for PDB %s, PMID: %s", pdbid, pmid)
            except IndexError, err:
                logging.warning("No results for PDB %s, PMID: %s", pdbid, pmid)
                issue    = ''
                date     = row['year']
                exc_l.append((pdbid, err))
            except KeyError, err:
                logging.warning("Missing PMC data for PDB %s, PMID: %s: %s", pdbid, pmid, str(err))
                issue    = ''
                date     = row['year']
                exc_l.append((pdbid, err))
            #else:
            #    logging.debug("Success retrieving PMC data for PDB %s, PMID: %s", pdbid, pmid)

            try:
                rowdict['source'] = "{journal_abbrev} {date}; {journal_volume}{issue}:{page_first}-{page_last} doi:{pdbx_database_id_doi}".format(
                    date=date, issue=issue, **row)
            except StandardError, err:
                logging.warning("Error constructing source field: %s\n%s", err, row.keys())

            citationwriter.writerow({k: v.replace('\t','  ').encode('utf8')
                                     if not isinstance(v,int) else v for k,v in rowdict.iteritems()})
            pdbs_completed += 1
            seen_pmid.add(int(pmid))

            if pdbs_completed % 5000 == 0:
                logging.info("Successfully retrieved citations for %d PDBIDs so far.", pdbs_completed)

        except KeyboardInterrupt:
            raise
        except Queue.Empty:
            logging.error("Timeout reached waiting for queue", exc_info=True)
            raise
        except StandardError, e:
            logging.error("Error processing PDB %s", pdbid, exc_info=True)
            exc_l.append((pdbid,e))
            pdbs_failed += 1
        else:
            sleep(0.005)
            pdbs_completed += 1

    logging.info("Dump worker finalising. PDBs written: %d, failed: %d, Duplicate PMIDs: %s",
                 pdbs_completed, pdbs_failed, pdbs_duped)
    citation_fh.close()
    return


def load_res_map(connection, resmapfields, timestamp=None):
    logging.info("Now loading res_map dump %s", RES_MAP)
    timestamp = timestamp or str(datetime.fromtimestamp(getmtime(RES_MAP)))
    with connection.begin() as transaction:
        connection.execute("COPY {}.res_map({}) FROM '{}'".format(SCHEMA, ','.join(resmapfields), RES_MAP))
        connection.execute("COMMENT ON TABLE {schema}.res_map IS 'Most recent data: {timestamp}'"\
                           .format(schema=SCHEMA, timestamp=str(timestamp)) )

def load_map_regions(connection, timestamp=None):
    fields = ','.join(['pdb','pdb_chain_id','sifts_res_num_start','sifts_res_num_end','db_source','db_accession_id','db_version','db_coord_sys'])
    logging.info("Now loading map_regions dump %s", MAP_REGIONS)
    timestamp = timestamp or datetime.fromtimestamp(getmtime(MAP_REGIONS))
    with connection.begin() as transaction:
        connection.execute("COPY {}.map_regions({}) FROM '{}'".format(SCHEMA, fields, MAP_REGIONS))
        connection.execute("COMMENT ON TABLE {schema}.map_regions IS 'Most recent data: {timestamp}'"\
                           .format(schema=SCHEMA, timestamp=str(timestamp)) )

def load_citations(connection, timestamp=None):
    logging.info("Now loading citations dump %s", CITATIONS)
    timestamp = timestamp or str(datetime.fromtimestamp(getmtime(CITATIONS)))
    with connection.begin() as transaction:
        connection.execute("COPY {}.citations({}) FROM '{}'".format(SCHEMA, ','.join(['pubmed_id','source','title','abstract']), CITATIONS))
        connection.execute("COMMENT ON TABLE {schema}.citations IS 'Most recent data: {timestamp}'"\
                           .format(schema=SCHEMA, timestamp=str(timestamp)) )
    logging.info("Finished loading citation dump")

    
def insert_ligands(connection, timestamp=None):
    logging.info("Inserting peptide ligands")
    try:
        connection.execute("""
                       WITH pss AS (
                           SELECT Structure_ID as structure_id, entity_id, pdb_strand_id as chain_id
                           FROM {mmcif}.pdbx_poly_seq_scheme
                           GROUP BY structure_id, entity_id, pdb_strand_id
                       )
                       INSERT INTO {schema}.peptide_ligands(pdb, pdb_chain_id, description, fragment, seq)
                       SELECT      e.structure_id as pdb, pss.chain_id as pdb_chain_id, e.pdbx_description as description, e.pdbx_fragment as fragment, ep.pdbx_seq_one_letter_code as seq
                       FROM        {mmcif}.entity e
                       LEFT JOIN   {mmcif}.entity_poly ep ON (e.structure_id = ep.structure_id AND e.id = ep.entity_id)
                       LEFT JOIN   pss ON (e.structure_id = pss.structure_id AND e.id = pss.entity_id)
                       WHERE       ep.type = 'polypeptide(L)' AND length(ep.pdbx_seq_one_letter_code_can) <= 10
                       ORDER BY    pdb, pdb_chain_id
                       """.format(schema=SCHEMA, mmcif=MMCIF))
    except SQLAlchemyError, err:
        logging.error("Error inserting peptide ligands:\n%s", str(err))
    else:
        if timestamp:
            connection.execute("COMMENT ON TABLE {schema}.peptide_ligands IS 'Most recent data: {timestamp}'"\
                               .format(schema=SCHEMA, timestamp=str(timestamp)) )

    logging.info("Inserting ligands")
    try:
        connection.execute("""
                       INSERT INTO {schema}.ligands(pdb, pdb_chain_id, res_num, ins_code, het_id, name)
                       SELECT      e.structure_id as pdb, nps.pdb_strand_id as pdb_chain_id, nps.pdb_seq_num::integer as res_num, nps.pdb_ins_code as ins_code, nps.pdb_mon_id as het_id,
                                   CASE WHEN enp."name" IS NOT NULL THEN enp.name ELSE e.pdbx_description END as name
                       FROM        {mmcif}.entity e
                       LEFT JOIN   {mmcif}.pdbx_entity_nonpoly enp ON (e.structure_id = enp.structure_id AND e.id = enp.entity_id)
                       JOIN        {mmcif}.pdbx_nonpoly_scheme nps ON (e.structure_id = nps.structure_id AND e.id = nps.entity_id)
                       WHERE       e.type != 'water'
                       ORDER BY    pdb, pdb_chain_id, res_num
                       """.format(schema=SCHEMA, mmcif=MMCIF))
    except SQLAlchemyError, err:
        logging.error("Error inserting ligands:\n%s", str(err))
    else:
        if timestamp:
            connection.execute("COMMENT ON TABLE {schema}.ligands IS 'Most recent data: {timestamp}'"\
                               .format(schema=SCHEMA, timestamp=str(timestamp)) )
        

def main(opt):
    '''
    '''
    resmapfields = ['entry','entity_id','sifts_res_num','sifts_res_name','pdb','pdb_chain_id','pdb_res_num',
                    'pdb_ins_code','pdb_res_name','uniprot','uniprot_res_num',
                    'uniprot_res_name','sstruct','is_observed','is_modified','is_conflicting']

    if all(isfile(f) for f in (RES_MAP, MAP_REGIONS, CITATIONS)) and opt.skip_dump and not opt.citations_only:
        logging.info("Dump file %s and %s already exist. Reusing.", RES_MAP, MAP_REGIONS)
        most_recent = datetime.fromtimestamp(getmtime(RES_MAP))
    else:
        # REGISTER NEW DIALECT TO WRITE TAB-DELIMITED FILES
        csv.register_dialect('tabs', delimiter='\t')

        if isfile(RES_MAP) or isfile(MAP_REGIONS):
            if opt.citations_only:
                resmapfh = open('/dev/null', 'w')
                mapregfh = open('/dev/null', 'w')
            elif opt.confirm or raw_input("Files exist. Overwrite? ").strip().lower() == 'y':
                resmapfh = open(RES_MAP, 'w')
                mapregfh = open(MAP_REGIONS, 'w')
            else:
                logging.warning("Not overwriting. Abort.")
                sys.exit(1)
        else:
            resmapfh = open(RES_MAP, 'w')
            mapregfh = open(MAP_REGIONS, 'w')

        # DICTIONARY WRITERS TO WRITE DATA TO FLAT FILES
        resmapwriter = csv.DictWriter(resmapfh, restval='\N', extrasaction='ignore',
                                      dialect='tabs',
                                      fieldnames = resmapfields)

        mapregwriter = csv.DictWriter(mapregfh, restval='\N', extrasaction='ignore',
                                      dialect='tabs',
                                      fieldnames=['pdb','pdb_chain_id','region_start','region_end','dbSource','dbAccessionId','dbVersion','dbCoordSys'])

        manager = mp.Manager()
        in_q    = mp.Queue()
        exc_l   = manager.list()
        dumproc = mp.Process(target=dump_citations, name="Citation-Process", args=(in_q, exc_l, connection))
        dumproc.start()

        FILES = os.listdir(SIFTS_DIR) if not opt.test else [xml for xml in os.listdir(SIFTS_DIR) if xml[:4].upper() in PDB_TEST_SET]

        # KEEP ONLY FILES FROM TEST SET
        # FILES = [xml for xml in FILES if xml[:4].upper() in PDB_TEST_SET]

        # INITIALIZE PROGRESSBAR
        bar = ProgressBar(widgets=['SIFTS: ', SimpleProgress(), ' ', Percentage(), Bar()], maxval=len(FILES)).start()

        # KEEP TRACK OF POSSIBLE RESIDUE ANNOTATIONS
        annotations = set() # [Conflict,Not_Observed]

        row_count = 0
        most_recent = datetime.fromtimestamp(0)
        for counter, filename in enumerate(FILES,1):
            bar.update(counter)

            if not filename.endswith('.gz'):
                continue

            pdb = filename.split('.',1)[0].upper()
            path = os.path.join(SIFTS_DIR, filename)

            # IGNORE EMPTY FILES
            if not isfile(path) or getsize(path) == 0:
                continue

            XML = gzip.open(path, 'rb')

            # PARSE FILE INTO XML TREE
            try:
                root = etree.parse(XML)
            except StandardError, err:
                logging.error("Error parsing compressed XML file %s: %s", path, str(err))
                continue
            else:
                in_q.put(pdb)
                #row_count += 1
                #if counter % 2500 == 0:
                #    logging.info("%8d entries processed. ", row_count)
                #continue
                if opt.citations_only:
                    continue

            entry = root.getroot()

            try:
                ns = '{%s}' % entry.nsmap[None]
                entities = entry.findall(ns+'entity')
            except AttributeError:
                entities = entry.findall('{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity')
                if not entities:
                    ns = '{http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd}'
                    entities = entry.findall(ns+'entity')
                    if not entities:
                        logging.warning("Found no XML entities on %s. Unknown namespace? ", filename)
                else:
                    ns = '{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}'
            finally:
                mapregionpath = '{0}listMapRegion/{0}mapRegion'.format(ns)
                
            dbversions = { dbelem.get('dbSource'): dbelem.get('dbVersion') for dbelem in entry.findall('{0}listDB/{0}db'.format(ns)) }
                
            # ITERATE THROUGH ALL ENTITIES (CHAINS)
            for entity in entities: #entry.findall(NS+'entity'):
                entity_id = entity.get('entityId')

                # ACTUAL PDB CHAIN ID CAN DIFFER FROM ENTITY ID / 1A4Y
                pdb_chain_id = None

                # ENTITY SEGMENTS
                for segment in entity.getchildren():

                    # LISTS OF RESIDUES
                    for listresidue in segment.findall(ns+'listResidue'): #NS+

                        # RESIDUES
                        for residue in listresidue.getchildren():

                            # SOME DEFAULT VALUES
                            row = {'sstruct': 'L', 'is_observed':1, 'is_modified':0, 'is_conflicting':0}

                            row['entry'] = pdb
                            row['entity_id'] = entity_id
                            row['sifts_res_num'] = residue.get('dbResNum')
                            row['sifts_res_name'] = residue.get('dbResName')

                            # THE ANNOTATIONS FOR EACH RESIDUE
                            for child in residue.getchildren():
                                if child.get('dbSource') == 'UniProt':
                                    row['uniprot'] = child.get('dbAccessionId').upper()
                                    row['uniprot_res_num'] = child.get('dbResNum')
                                    row['uniprot_res_name'] = child.get('dbResName')

                                elif child.get('dbSource') == 'PDB':
                                    row['pdb'] = child.get('dbAccessionId').upper()
                                    row['pdb_chain_id'] = child.get('dbChainId')
                                    row['pdb_res_name'] = child.get('dbResName')

                                    pdb_chain_id = row['pdb_chain_id']

                                    # GET THE INSERTION CODE (CONCATENATED IN SIFTS)
                                    dbresnum = child.get('dbResNum')

                                    if dbresnum[-1].isdigit():
                                        row['pdb_res_num'] = dbresnum
                                        row['pdb_ins_code'] = ' '

                                    else:
                                        row['pdb_res_num'] = dbresnum[:-1]
                                        row['pdb_ins_code'] = dbresnum[-1]

                                # GET THE SECONDARY STRUCTURE CODE
                                elif child.get('property') == 'codeSecondaryStructure':
                                    row['sstruct'] = child.text

                                # EXTRACT OTHER MSD RESIDUE ANNOTATION
                                elif child.get('property') == 'Annotation':
                                    annotation = child.text.replace('\n          ',' ').strip()

                                    if annotation == 'Not_Observed': row['is_observed'] = '0'
                                    elif annotation == 'PDB modified': row['is_modified'] = '1'
                                    elif annotation == 'Conflict': row['is_conflicting'] = '1'
                                    else: annotations.add(annotation)

                            # WRITE RESIDUE INFORMATION TO FILE
                            resmapwriter.writerow(row)
                            row_count += 1

                    # GET ALL THE MAPPED REGIONS OF THE SEGMENT
                    for mapregion in segment.findall(mapregionpath):
                        start, end = mapregion.get('start'), mapregion.get('end')

                        # IGNORE CASES WHERE THERE IS NO DEFINED REGION
                        if not start and not end:
                            continue

                        # ITERATE THROUGH THE EXTERNAL DATABASE MAPPINGS
                        for db in mapregion.iter(ns+'db'):
                            row = db.attrib

                            # USE THE ORIGINAL PDB CHAIN ID HERE
                            row.update({'pdb': pdb, 'pdb_chain_id': pdb_chain_id, 'region_start': start,'region_end': end})
                            
                            # ADD DB VERSIONS
                            row['dbVersion'] = dbversions.get(row['dbSource'], '\N')

                            mapregwriter.writerow(row)

            xml_mtime = datetime.fromtimestamp(getmtime(path))
            if xml_mtime > most_recent:
                most_recent = xml_mtime

            if counter % 2500 == 0:
                logging.info("%8d rows written.", row_count)
        else:
            bar.finish()
            resmapfh.close()
            mapregfh.close()
            logging.info("Finished processing SIFTS. Waiting for citations process...")

            in_q.put(None)
            in_q.close()
            try:
                dumproc.join(12000 if opt.citations_only else 3600)
            except mp.TimeoutError:
                logging.warn("Timeout processing SIFTS. Continuing...")
            else:
                logging.info("Finished processing SIFTS.")


        print annotations

    if not opt.skip_load and (opt.confirm or raw_input(("Tables in schema '{}' are about be dropped and recreated. " +
                                                        "Are you sure you want to continue? [y/N] ").format(SCHEMA)).strip().lower() == 'y'):
        logging.info("Creating tables...")
        res_map, map_regions, citations = create_tables(metadata)
    elif opt.skip_load:
        res_map, map_regions, citations = reflect_tables(metadata)
    else:
        print >> sys.stderr, "Aborting."
        sys.exit(1)


    ## RES_MAP dump loading and index creation
    if not opt.skip_load:
        try:
            load_res_map(connection, resmapfields, timestamp=most_recent)
        except SQLAlchemyError, err:
            logging.error("Error loading res_map table!\n%s", str(err))
    try:
        Index('idx_res_map_entry', res_map.c.entry, res_map.c.entity_id, res_map.c.sifts_res_num).create(engine) # NOT UNIQUE BECAUSE A MODIFIED RESIDUE CAN SPAN MORE THAN ONE RESIDUE -> 1GVW A 54
        Index('idx_res_map_pdb', res_map.c.pdb, res_map.c.pdb_chain_id, res_map.c.pdb_res_num).create(engine)
        Index('idx_res_map_uniprot', res_map.c.uniprot, res_map.c.uniprot_res_num).create(engine)
        Index('idx_res_map_sstruct', res_map.c.pdb, res_map.c.pdb_chain_id, res_map.c.sstruct_serial).create(engine)
    except ProgrammingError:
        logging.warning("Indices on res_map already exist.")


    ## MAP_REGIONS dump loading and index creation
    if not opt.skip_load:
        try:
            load_map_regions(connection,  timestamp=most_recent)
        except SQLAlchemyError, err:
            logging.error("Error loading map_regions:\n%s", str(err))
    try:
        Index('idx_map_regions_pdb', map_regions.c.pdb, map_regions.c.pdb_chain_id).create(engine)
        Index('idx_map_regions_dbref', map_regions.c.db_source, map_regions.c.db_accession_id).create(engine)
    except ProgrammingError:
        logging.warning("Indices on map_regions already exist.")


    ## CITATIONS dump loading and index creation
    if not opt.skip_load:
        try:
            load_citations(connection, timestamp=most_recent)
        except SQLAlchemyError, err:
            logging.error("Error loading citations:\n%s", str(err))
    try:
        Index('idx_citations_abstract_tsvector', func.to_tsvector('english', citations.c.abstract), postgresql_using="gin").create(engine)
    except ProgrammingError, err:
        logging.warning("Index on citations already exists: %s", err)


    update_sstruct_serials(connection)
    insert_disordered_regions(connection)
    insert_prot_fragments(connection, timestamp=most_recent)
    insert_ligands(connection, timestamp=most_recent)

    connection.close()

if __name__ == '__main__':
    main(opt)
