#!/usr/bin/python
import os
import csv
import json
import warnings

from sqlalchemy import *
from sqlalchemy.event import listen
from sqlalchemy.engine.url import URL
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemyutils import Cube, ArrayXd

from eyesopen.oechem import reset_charges, OEDots, oemolistream, OECount, OEIsHeavy, OEGetSDData, OESuppressHydrogens
from eyesopen.oemolprop import nprs
from usrcat.toolkits.oe import generate_moments

# REGISTER NEW DIALECT TO WRITE TAB-DELIMITED FILES
csv.register_dialect('tabs', delimiter='\t')

PDBCHEM_CONF_DIR    = '/tlbnas/store/pdbe/pdbechem/conformers'
CHEM_CONFORMER_DUMP = '/tlbnas/temp/bahamut/chem_comp_conformers.pdbchem'

engine      = create_engine(URL(drivername='postgresql+psycopg2', username='adrian',
                                password='1r1d1Um', host='bahamut.bioc.cam.ac.uk',
                                port=5432, database='cryst'),
                            execution_options={'autocommit':True}, echo=False)
connection  = engine.connect()
metadata    = MetaData(bind=connection)

def create_table():
    '''
    '''
    table = Table('chem_comp_conformers', metadata,
                  Column('het_id', String(3), nullable=False, primary_key=True),
                  Column('conformer', Integer, nullable=False, autoincrement=False, primary_key=True),
                  Column('energy', Float(8,2), nullable=False),
                  Column('npr1', Float(4,3)),
                  Column('npr2', Float(4,3)),
                  Column('usr_space', Cube, nullable=False),
                  Column('usr_moments', ArrayXd, nullable=False),
                  schema='pdbchem')

    # MAKE THESE TABLES PUBLICLY ACCESSIBLE
    listen(table, "after_create", DDL("GRANT SELECT ON TABLE pdbchem.chem_comp_conformers TO public"))

    table.drop(checkfirst=True)
    table.create(checkfirst=True)

def main():
    '''
    '''
    ofs = open(CHEM_CONFORMER_DUMP,'w')
    writer = csv.writer(ofs, dialect='tabs')

    dots = OEDots(100,1, "Chemical components")

    ifs = oemolistream()

    for file in os.listdir(PDBCHEM_CONF_DIR):
        path = os.path.join(PDBCHEM_CONF_DIR, file)

        if os.path.getsize(path) <= 15: continue

        ifs.open(str(path))

        try:
            molecule = ifs.GetOEMols().next()
        except StopIteration:
            print file
            continue

        # do not generate moments for small structures
        if OECount(molecule, OEIsHeavy()) < 7: continue

        # to be compatible with structures from the pdb
        OESuppressHydrogens(molecule)

        # get the USRCAT moments for all conformers of the molecule in one go
        moments = generate_moments(molecule)

        for conformer in molecule.GetConfs():
            title = conformer.GetTitle()

            if title.find('_') == -1: continue

            het_id, conf_id = title.rsplit('_', 1)
            energy = OEGetSDData(conformer, 'mmff94s_NoEstat')

            confmoments = moments[conformer.GetIdx()]

            usrspace = '(' + ','.join(map(str,confmoments[:12])) + ')'
            usrmoments = '{' + ','.join(map(str,confmoments)) + '}'

            try:
                npr1, npr2 = nprs(conformer)
            except (AssertionError, ValueError) as e:
                npr1 = npr2 = None

            writer.writerow([het_id, conf_id, energy,  npr1, npr2, usrspace, usrmoments])

        dots.Update()

    dots.Total()
    ofs.close()

    create_table()

    connection.execute("COPY pdbchem.chem_comp_conformers FROM '{}'".format(CHEM_CONFORMER_DUMP))

    DDL("CREATE INDEX idx_chem_comp_conformers_usr_space ON pdbchem.chem_comp_conformers USING GIST(usr_space) WITH (FILLFACTOR=100)").execute(bind=engine)
    DDL("ALTER TABLE pdbchem.chem_comp_conformers CLUSTER ON idx_chem_comp_conformers_usr_space").execute(bind=engine)

    # cluster table using the gist index
    connection.execute("CLUSTER pdbchem.chem_comp_conformers")

main()
