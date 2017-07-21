#!/usr/bin/env python
import os, sys
import re
import logging

from datetime import datetime
from lxml import etree

from sqlalchemy                 import create_engine, Column, Table, Index, DDL
from sqlalchemy.engine.url      import URL
from sqlalchemy.orm             import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.types           import Integer, String, Text,  TypeDecorator
from sqlalchemy.sql.expression  import cast, func
from sqlalchemy.exc             import *


def user_confirm(msg='Irreversible stuff will happen!'):
    return raw_input(msg + " Are you sure you want to do this? [y/N] ").strip().lower() == 'y'

class OBUniversalSmiles(TypeDecorator):
    impl = Text

    def bind_expression(self, value):
        try:
            return func.openbabel.convert_molecule(value, 'SMI', 'SMI', 'U')
        except Exception, err: # Not sure an exception is raised at this point. Just in case...
            logging.warning("%sFailed to standardise SMILES '%s': %s%s", Colors.WARNING, value, err, Colors.ENDC)
            return value
        else:
            logging.debug("Successfully converted %s to OpenBabel Universal SMILES", value)
       
            
class OECanonicalSmiles(TypeDecorator):
    impl = Text

    def bind_expression(self, value):
        try:
            return func.openeye.standardise_smiles(value)
        except Exception, err: # Not sure an exception is raised at this point. Just in case...
            logging.warning("%sFailed to standardise SMILES '%s': %s%s", Colors.WARNING, value, err, Colors.ENDC)
            return value
        else:
            logging.debug("Successfully converted %s to OpenEye canonical SMILES", value)


## Command-line & database initialization
if __name__ == '__main__':
    import argparse
    from getpass  import getpass

    parser = argparse.ArgumentParser()
    parser.add_argument('emolecules', nargs='?', type=argparse.FileType('r'), help="Emolecules SMI file",
                        default='emolecules.smi' if os.path.exists('emolecules.smi') else None)
    parser.add_argument('-s', '--start-line', dest="start_line", type=int, help="Line to start processing at.", default=None)
    parser.add_argument('-u', '--user',   dest="dbuser", help="Username on database", default='bernardo')
    parser.add_argument('-H', '--host',   dest="dbhost", help="Location of the database server/file", default='192.168.1.86')
    parser.add_argument('-p', '--port',   dest="dbport", type=int, help="Port of the database server", default=5432)
    parser.add_argument('-D', '--db',     dest="dbname", help="Name of the database", default='cryst')
    parser.add_argument('-S', '--schema', dest="schema", help="Schema to use on the database", default='emolecules')
    parser.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Echo SQL commands", default=False)
    parser.add_argument('-L', '--loglevel',dest="loglevel", action='store', choices=['debug','info','warn','error'],
                        help="Logging level to display (default: information)", default='warn')
    
    opt = parser.parse_args()
    optdict = vars(opt)
    logging.basicConfig(level=getattr(logging, optdict.get('loglevel').upper()))

    engine = create_engine(URL(drivername='postgresql+psycopg2', username=opt.dbuser, 
                               password=getpass("DB Password for '%s' on %s:" % (opt.dbuser, opt.dbhost)),
                               host=opt.dbhost, port=opt.dbport, database=opt.dbname),
                           echo_pool=opt.verbose, pool_size=10, pool_recycle=600)
    Base = declarative_base(bind=engine)


class EMolecule(Base):
    __table__ = Table('molecules', Base.metadata,
                      Column('ism_emol', Text, primary_key=True),
                      Column('ism_ob',   OBUniversalSmiles, index=True),
                      Column('ism_oe',   OECanonicalSmiles, index=True),
                      Column('version_id', Integer, nullable=False),
                      Column('parent_id',  Integer, nullable=False),
                      schema=opt.schema
                      )

    def __init__(self, smiles, version_id, parent_id):
        self.ism_emol = smiles
        self.ism_ob   = smiles
        self.ism_oe   = smiles
        self.version_id = version_id
        self.parent_id  = parent_id
        

def main(smiles_file, start_line=None):
    if user_confirm("Tables will be dropped! "):
        logging.info("Dropping tables on schema %s...", opt.schema)
        Base.metadata.drop_all(checkfirst=True)

    logging.info("Creating tables on schema %s...", opt.schema)
    Base.metadata.create_all(checkfirst=True)

    Session  = sessionmaker(bind=engine)
    session  = Session()

    logging.info("Connected to database. Now parsing...")
    
    for i, line in enumerate(smiles_file, start=1):
        try:
            if start_line:
                if i < start_line:
                    continue
                elif i == start_line:
                    logging.info("Starting processing at line %d", i)
            
            smiles, version_id, parent_id = line.split()
            if version_id == 'version_id':
                continue
            emol = session.merge(EMolecule(smiles, version_id, parent_id))
        
            if i % 10000 == 0:
                logging.info("%d molecules processed. Committing.", i)
                session.commit()
            elif i % 1000 == 0:
                logging.info("%d molecules processed", i)
                session.flush()
        except (SQLAlchemyError, StandardError), err:
            logging.error(err.message.strip())
            session.rollback()
    else:
        session.commit()
        logging.info("%d molecules loaded in total" % i)

    
    try:
        mtime = os.path.getmtime(smiles_file.name)
        timestamp = str(datetime.fromtimestamp(mtime))
        session.execute("COMMENT ON TABLE {schema}.molecules IS 'Last modification time: {timestamp}'"\
                           .format(schema=opt.schema, timestamp=timestamp) )
    except AttributeError:
        logging.warning("Data provided by pipe, file modification time not available")
        
    session.commit()


if __name__ == '__main__':
    main(opt.emolecules, opt.start_line)
