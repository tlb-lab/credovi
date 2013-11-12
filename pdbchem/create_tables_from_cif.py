import os
import json
from subprocess import Popen, PIPE

from sqlalchemy import Column, Table, Integer, String, create_engine, MetaData, Date, Text, Numeric, Index
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY
from Bio.PDB import MMCIF2Dict
from progressbar import ProgressBar, Percentage, Bar, SimpleProgress

LIG_CIF_DIR         = '/tlbnas/mirror/pdbe/pdbechem/files/mmcif'
PDB_CHEM_COMPS_DUMP = '/tlbnas/temp/bahamut/pdb_chem_comps.pdbchem'

engine      = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                port=5432, database='cryst'),
                            execution_options={'autocommit':True}, echo=False)
connection  = engine.connect()
metadata    = MetaData(bind=connection)

def create_tables():
    '''
    '''
    pdb_chem_comps = Table('pdb_chem_comps', metadata,
                           Column('id', String(3), nullable=False, primary_key=True),
                           Column('name', Text),
                           Column('monomer_type', String(32)),
                           Column('pdbx_type', String(6)),
                           Column('formula', String(48)),
                           Column('mon_nstd_parent_comp_id', Text),
                           Column('pdbx_inchi', Text),
                           Column('pdbx_inchikey', Text),
                           Column('pdbx_synonyms', Text),
                           Column('pdbx_formal_charge', Integer),
                           Column('pdbx_initial_date', Date),
                           Column('pdbx_modified_date', Date),
                           Column('pdbx_ambiguous_flag', String(1)),
                           Column('pdbx_release_status', String(8)),
                           Column('pdbx_replaced_by', String(3)),
                           Column('pdbx_replaces', String(13)),
                           Column('formula_weight', Numeric(8,3)),
                           Column('one_letter_code', String(4)),
                           Column('three_letter_code', String(3)),
                           Column('pdbx_model_coordinates_details', Text),
                           Column('pdbx_model_coordinates_missing_flag', String(1)),
                           Column('pdbx_ideal_coordinates_details', Text),
                           Column('pdbx_ideal_coordinates_missing_flag', String(1)),
                           Column('pdbx_model_coordinates_db_code', String(19)),
                           Column('pdbx_subcomponent_list', ARRAY(String(3))),
                           Column('pdbx_processing_site', String(4)),
                           schema='pdbchem')

    Index('idx_pdb_chem_comps_mon_nstd_parent_comp_id', pdb_chem_comps.c.mon_nstd_parent_comp_id)

    non_std_res = Table('non_std_res', metadata,
                        Column('het_id', String(3), nullable=False, primary_key=True),
                        Column('std_het_id', String(3), nullable=False, primary_key=True),
                        schema='pdbchem')

    Index('idx_non_std_res_std_het_id', non_std_res.c.std_het_id, non_std_res.c.het_id, unique=True)

    metadata.drop_all(checkfirst=True)
    metadata.create_all(checkfirst=True)

def main():
    '''
    '''
    out = open(PDB_CHEM_COMPS_DUMP,'w')

    cifs = os.listdir(LIG_CIF_DIR)
    cifs.sort()

    # INITIALIZE PROGRESSBAR
    bar = ProgressBar(widgets=['CIF files: ', SimpleProgress(), ' ', Percentage(), Bar()], maxval=len(cifs)).start()

    for counter, cif in enumerate(cifs,1):
        bar.update(counter)

        path = os.path.join(LIG_CIF_DIR,cif)

        if path[-3:] != 'cif': continue

        data = MMCIF2Dict.MMCIF2Dict(path)

        fields = ['\N'] * 26

        fields[0] = data.get('_chem_comp.id', '\N')
        fields[1] = data.get('_chem_comp.name', '\N')
        fields[2] = data.get('_chem_comp.type', '\N')
        fields[3] = data.get('_chem_comp.pdbx_type', '\N')
        fields[4] = data.get('_chem_comp.formula', '\N')
        fields[5] = data.get('_chem_comp.mon_nstd_parent_comp_id', '\N')

        #fields[8] = data.get('_chem_comp.pdbx_synonyms', '\N')
        fields[9] = data.get('_chem_comp.pdbx_formal_charge', '\N')
        fields[10] = data.get('_chem_comp.pdbx_initial_date', '\N')
        fields[11] = data.get('_chem_comp.pdbx_modified_date', '\N')
        fields[12] = data.get('_chem_comp.pdbx_ambiguous_flag', '\N')
        fields[13] = data.get('_chem_comp.pdbx_release_status', '\N')
        fields[14] = data.get('_chem_comp.pdbx_replaced_by', '\N')
        fields[15] = data.get('_chem_comp.pdbx_replaces', '\N')
        fields[16] = data.get('_chem_comp.formula_weight', '\N')
        fields[17] = data.get('_chem_comp.one_letter_code', '\N')
        fields[18] = data.get('_chem_comp.three_letter_code', '\N')
        fields[19] = data.get('_chem_comp.pdbx_model_coordinates_details', '\N')
        fields[20] = data.get('_chem_comp.pdbx_model_coordinates_missing_flag', '\N')
        fields[21] = data.get('_chem_comp.pdbx_ideal_coordinates_details', '\N')
        fields[22] = data.get('_chem_comp.pdbx_ideal_coordinates_missing_flag', '\N')
        fields[23] = data.get('_chem_comp.pdbx_model_coordinates_db_code', '\N')

        subcomponents = data.get('_chem_comp.pdbx_subcomponent_list', '\N')

        if subcomponents != '?':
            subcomponents = subcomponents.strip('"').split(' ')
            subcomponents = '{' + ','.join(subcomponents) + '}'

        fields[24] = subcomponents
        fields[25] = data.get('_chem_comp.pdbx_processing_site', '\N')

        try:
            fields[6] = data['_pdbx_chem_comp_descriptor.descriptor'][data['_pdbx_chem_comp_descriptor.type'].index('InChI')]
            fields[7] = data['_pdbx_chem_comp_descriptor.descriptor'][data['_pdbx_chem_comp_descriptor.type'].index('InChIKey')]
        except:
            pass

        # REMOVE NEWLINES AND REPLACE QUESTION MARKS
        fields = [field.strip().replace('\n','').replace('?','\N') for field in fields]
        fields = [field.replace('"','') for field in fields]

        line = '\t'.join(fields)

        print >> out, line

    bar.finish()
    out.close()

    create_tables()

    connection.execute("COPY pdbchem.pdb_chem_comps FROM '{}'".format(PDB_CHEM_COMPS_DUMP))

    connection.execute("""
                        INSERT  INTO pdbchem.non_std_res
                        SELECT  sq.het_id, sq.std_het_id
                        FROM    (
                                SELECT       DISTINCT id as het_id, unnest(string_to_array(mon_nstd_parent_comp_id, ',')) as std_het_id
                                FROM         pdbchem.pdb_chem_comps
                                WHERE        mon_nstd_parent_comp_id IS NOT NULL
                                ORDER BY     1,2
                                ) sq
                        WHERE   LENGTH(sq.std_het_id) <= 3
                       """)

main()
