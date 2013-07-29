import os
import gzip
import sys
import json
import csv
from copy import copy
import xml.etree.cElementTree as ctree

from sqlalchemy import *
from sqlalchemy.event import listen
from sqlalchemy.schema import CreateIndex
from sqlalchemy.engine.url import URL
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

engine      = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                port=5432, database='cryst'),
                            execution_options={'autocommit':True}, echo=False)
connection  = engine.connect()
metadata    = MetaData(bind=connection)

CONFIG      = json.loads(open('/home/adrian/Dropbox/Research/Configurations/project.json').read())
CREDO_CONF  = json.loads(open('/home/adrian/Dropbox/Research/Configurations/credo.json').read())

SIFTS_DIR   = '/tlbnas/mirror/pdbe/sifts/xml'
PDB_TEST_SET= frozenset(CREDO_CONF['PDB TEST SET'])

SCHEMA      = 'pdb'
NS         = '{http://www.efamily.org.uk/xml/efamily/2004/08/14/eFamily.xsd}'

RES_MAP  = '/tlbnas/temp/bahamut/res_map.pdb'
MAP_REGIONS = '/tlbnas/temp/bahamut/map_regions.pdb'

MAPREGIONPATH = '{0}listMapRegion/{0}mapRegion'.format(NS)

def create_tables():
    '''
    '''
    res_map = Table('res_map', metadata,
                    Column('res_map_id', Integer, primary_key=True),
                    Column('entry', String(4), nullable=False),
                    Column('entity_id', String(2), nullable=False),
                    Column('sifts_res_num', Integer, nullable=False),
                    Column('sifts_res_name', String(3)),
                    Column('pdb', String(4)),
                    Column('pdb_chain_id', String(1)),
                    Column('pdb_res_num', Integer),
                    Column('pdb_ins_code', String(1)),
                    Column('pdb_res_name', String(3)),
                    Column('uniprot', String(6)),
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
                        Column('pdb_chain_id', String(1), nullable=False),
                        Column('sifts_res_num_start', Integer, nullable=False),
                        Column('sifts_res_num_end', Integer, nullable=False),
                        Column('db_source', String(12)),
                        Column('db_accession_id', String(12)),
                        Column('db_version', String(12)),
                        Column('db_coord_sys', String(12)),
                        schema=SCHEMA)

    disordered_regions = Table('disordered_regions', metadata,
                               Column('map_region_id', Integer, primary_key=True),
                               Column('pdb', String(4), nullable=False),
                               Column('pdb_chain_id', String(1), nullable=False),
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
                               Column('pdb_chain_id', String(1), nullable=False),
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

    uniprot_residues = Table('uniprot_residues', metadata,
                             Column('uniprot_residue_id', Integer, primary_key=True),
                             Column('uniprot', String(6)),
                             Column('uniprot_res_num', Integer),
                             Column('uniprot_res_name', String(1)),
                             schema=SCHEMA)

    # MAKE THESE TABLES PUBLICLY ACCESSIBLE
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE pdb.res_map TO public"))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE pdb.map_regions TO public"))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE pdb.disordered_regions TO public"))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE pdb.pdb_prot_fragments TO public"))
    listen(metadata, "after_create", DDL("GRANT SELECT ON TABLE pdb.pdb_prot_fragment_to_residue TO public"))

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

    return metadata

def update_sstruct_serials():
    '''
    '''
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM pdb.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        UPDATE  pdb.res_map m
                        SET     sstruct_serial = sq.sstruct_serial
                        FROM    (
                                WITH T1 AS
                                (
                                    SELECT  res_map_id, entry, entity_id, sifts_res_num, sstruct,
                                            LAG(sstruct,1,sstruct)
                                                OVER (ORDER BY entry, entity_id, sifts_res_num)
                                                IS DISTINCT FROM sstruct AS changes
                                    FROM pdb.res_map WHERE pdb = $1
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
                """

    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

def insert_disordered_regions():
    '''
    '''
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM pdb.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        INSERT INTO pdb.disordered_regions(pdb, pdb_chain_id, disordered_region_serial, region_seq, region_length, region_start, region_end, is_nterminus, is_cterminus)
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
                                    FROM    pdb.res_map
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
                """

    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

def insert_prot_fragments():
    '''
    '''
    # CREATE UNIQUE REPRESENTATION FOR PROTEIN FRAGMENTS IN THE PDB
    statement = """
                DO $$
                DECLARE
                    pdbcode character varying(4);
                BEGIN
                -- LOOP THROUGH ALL PDB ENTRIES
                    FOR pdbcode IN SELECT DISTINCT pdb FROM pdb.res_map WHERE pdb IS NOT NULL LOOP
                        EXECUTE
                        '
                        INSERT      INTO pdb.pdb_prot_fragments(pdb, pdb_chain_id, sstruct_serial, sstruct, fragment_size, fragment_seq, fragment_nterm, fragment_cterm)
                        SELECT      pdb, pdb_chain_id, sstruct_serial, sstruct,
                                    LENGTH(ARRAY_TO_STRING(ARRAY_AGG(bio.to_one_letter_code(sifts_res_name)),'''')) AS fragment_size,
                                    ARRAY_TO_STRING(array_agg(bio.to_one_letter_code(sifts_res_name) ORDER BY pdb_res_num),'''') AS seq,
                                    LAG(m.sstruct_serial) OVER(PARTITION BY pdb) AS nterminal,
                                    LEAD(m.sstruct_serial) OVER(PARTITION BY pdb) AS cterminal
                        FROM        pdb.res_map m
                        WHERE       m.pdb = $1
                        GROUP BY    pdb, pdb_chain_id, sstruct_serial, sstruct
                        ORDER BY    pdb, pdb_chain_id, sstruct_serial
                        ' USING pdbcode;

                    END LOOP;
                END$$;
                """

    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

    # CREATE A MAPPING BETWEEN PDB PROTEIN FRAGMENTS AND THE RESIDUES FROM RESMAP
    statement = """
                INSERT      INTO pdb.pdb_prot_fragment_to_residue
                SELECT      f.pdb_prot_fragment_id, m.res_map_id
                FROM        pdb.pdb_prot_fragments f
                JOIN        pdb.res_map m
                            ON m.pdb = f.pdb
                            AND m.pdb_chain_id = f.pdb_chain_id
                            AND m.sstruct_serial = f.sstruct_serial
                ORDER BY    f.pdb_prot_fragment_id, m.res_map_id
                """

    transaction = connection.begin()
    connection.execute(statement)
    transaction.commit()

def main():
    '''
    '''
    resmapfields = ['entry','entity_id','sifts_res_num','sifts_res_name','pdb','pdb_chain_id','pdb_res_num',
                    'pdb_ins_code','pdb_res_name','uniprot','uniprot_res_num',
                    'uniprot_res_name','sstruct','is_observed','is_modified','is_conflicting']

    # REGISTER NEW DIALECT TO WRITE TAB-DELIMITED FILES
    csv.register_dialect('tabs', delimiter='\t')

    resmapfh = open(RES_MAP, 'w')
    mapregfh = open(MAP_REGIONS, 'w')

    # DICTIONARY WRITERS TO WRITE DATA TO FLAT FILES
    resmapwriter = csv.DictWriter(resmapfh, restval='\N', extrasaction='ignore',
                                  dialect='tabs',
                                  fieldnames = resmapfields)

    mapregwriter = csv.DictWriter(mapregfh, restval='\N', extrasaction='ignore',
                                  dialect='tabs',
                                  fieldnames=['pdb','pdb_chain_id','region_start','region_end','dbSource','dbAccessionId','dbVersion','dbCoordSys'])

    FILES = os.listdir(SIFTS_DIR)

    # KEEP ONLY FILES FROM TEST SET
    # FILES = [xml for xml in FILES if xml[:4].upper() in PDB_TEST_SET]

    # INITIALIZE PROGRESSBAR
    bar = ProgressBar(widgets=['SIFTS: ', SimpleProgress(), ' ', Percentage(), Bar()], maxval=len(FILES)).start()

    # KEEP TRACK OF POSSIBLE RESIDUE ANNOTATIONS
    annotations = set() # [Conflict,Not_Observed]

    for counter, file in enumerate(FILES,1):
        bar.update(counter)

        if not file.endswith('.gz'): continue

        pdb = file.split('.',1)[0].upper()
        path = os.path.join(SIFTS_DIR, file)

        # IGNORE EMPTY FILES
        if os.path.getsize(path) == 0: continue

        XML = gzip.open(path, 'rb')

        # PARSE FILE INTO XML TREE
        root = ctree.parse(XML)
        entry = root.getroot()

        # ITERATE THROUGH ALL ENTITIES (CHAINS)
        for entity in entry.findall(NS+'entity'):
            entity_id = entity.get('entityId')

            # ACTUAL PDB CHAIN ID CAN DIFFER FROM ENTITY ID / 1A4Y
            pdb_chain_id = None

            # ENTITY SEGMENTS
            for segment in entity.getchildren():

                # LISTS OF RESIDUES
                for listresidue in segment.findall(NS+'listResidue'):

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

                # GET ALL THE MAPPED REGIONS OF THE SEGMENT
                for mapregion in segment.findall(MAPREGIONPATH):
                    start, end = mapregion.get('start'), mapregion.get('end')

                    # IGNORE CASES WHERE THERE IS NO DEFINED REGION
                    if not start and not end: continue

                    # ITERATE THROUGH THE EXTERNAL DATABASE MAPPINGS
                    for db in mapregion.iter(NS+'db'):
                        row = db.attrib

                        # USE THE ORIGINAL PDB CHAIN ID HERE
                        row.update({'pdb': pdb, 'pdb_chain_id': pdb_chain_id, 'region_start': start,'region_end': end})

                        mapregwriter.writerow(row)

    bar.finish()

    print annotations

    resmapfh.close()
    mapregfh.close()

    create_tables()

    res_map = metadata.tables['pdb.res_map']
    map_regions = metadata.tables['pdb.map_regions']

    transaction = connection.begin()
    connection.execute("COPY {}.res_map({}) FROM '{}'".format(SCHEMA, ','.join(resmapfields), RES_MAP))
    transaction.commit()

    # CREATE INDEXES ON RES MAP
    Index('idx_res_map_entry', res_map.c.entry, res_map.c.entity_id, res_map.c.sifts_res_num).create(engine) # NOT UNIQUE BECAUSE A MODIFIED RESIDUE CAN SPAN MORE THAN ONE RESIDUE -> 1GVW A 54
    Index('idx_res_map_pdb', res_map.c.pdb, res_map.c.pdb_chain_id, res_map.c.pdb_res_num).create(engine)
    Index('idx_res_map_uniprot', res_map.c.uniprot, res_map.c.uniprot_res_num).create(engine)
    Index('idx_res_map_sstruct', res_map.c.pdb, res_map.c.pdb_chain_id, res_map.c.sstruct_serial).create(engine)

    fields = ','.join(['pdb','pdb_chain_id','sifts_res_num_start','sifts_res_num_end','db_source','db_accession_id','db_version','db_coord_sys'])
    transaction = connection.begin()
    connection.execute("COPY {}.map_regions({}) FROM '{}'".format(SCHEMA, fields, MAP_REGIONS))
    transaction.commit()

    # CREATE INDEXES ON MAP REGIONS
    Index('idx_map_regions_pdb', map_regions.c.pdb, map_regions.c.pdb_chain_id).create(engine)
    Index('idx_map_regions_dbref', map_regions.c.db_source, map_regions.c.db_accession_id).create(engine)

    update_sstruct_serials()
    insert_disordered_regions()

    insert_prot_fragments()

    connection.close()

main()
