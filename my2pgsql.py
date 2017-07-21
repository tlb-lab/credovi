import os
import re
import sys
import logging
from getpass import getpass, getuser
from subprocess import Popen, PIPE

from sqlalchemy import *
from sqlalchemy.dialects.postgresql import TIMESTAMP
from sqlalchemy.dialects import mysql

# LOGGING
logger = logging.getLogger("my2pgsql.py")
logger.setLevel(logging.INFO)

# CREATE CONSOLE HANDLER AND SET LEVEL TO DEBUG
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# CREATE FORMATTER
formatter = logging.Formatter("%(levelname)s - %(message)s")

# ADD FORMATTER TO CH
ch.setFormatter(formatter)

# ADD CH TO LOGGER
logger.addHandler(ch)

def convert(name):
    '''
    Convert CamelCase to under_score. From StackOverflow.
    '''
    s = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s).lower()

def parse_options():
    '''
    '''
    from getpass  import getpass
    from optparse import OptionParser
    # PARSE COMMAND LINE
    usage  = "%prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-D", "--debug",
                      action  = "store_true",
                      dest    = "DEBUG",
                      default = False,
                      help    = 'Set logging level to debug and print more verbose output.')

    parser.add_option("--mysql-user",
                      dest    = "mysqluser",
                      default = getuser(), # None,
                      help    = "The MySQL user account to be used.")

    parser.add_option("--mysql-passwd",
                      dest    = "mysqlpasswd",
                      default = '',
                      help    = "The password for the given MySQL user account.")

    parser.add_option("--mysql-host",
                      dest    = "mysqlhost",
                      default = 'localhost',
                      help    = "The host of the MySQL database.")

    parser.add_option("--mysql-port",
                      dest    = "mysqlport",
                      default = 3306,
                      help    = "The port of the MySQL database.")

    parser.add_option("--mysql-db",
                      dest    = "mysqldb",
                      default = None,
                      help    = "The MySQL database to be migrated.")

    parser.add_option("--pgsql-user",
                      dest    = "pguser",
                      default = getuser(), #None,
                      help    = "The PosgreSQL user account to be used.")

    parser.add_option("--pgsql-passwd",
                      dest    = "pgpasswd",
                      default = None,
                      help    = "The password for the given PostgreSQL user account.")

    parser.add_option("--pgsql-host",
                      dest    = "pghost",
                      default = 'localhost',
                      help    = "The host of the PostgreSQL database.")

    parser.add_option("--pgsql-port",
                      dest    = "pgport",
                      default = 5432,
                      help    = "The port of the PostgreSQL database.")

    parser.add_option("--pgsql-db",
                      dest    = "pgdb",
                      default = None,
                      help    = "The target PostgreSQL database.")

    parser.add_option("--pgsql-schema",
                      dest    = "pgschema",
                      default = None,
                      help    = "The schema for the target PostgreSQL.")

    parser.add_option("--data",
                      action  = "store_true",
                      dest    = "data",
                      default = False,
                      help    = 'Migrate data.')

    parser.add_option("--tables",
                      dest    = "tables",
                      default = None,
                      help    = "The comma-separated list of MySQL tables that should be migrated.")

    parser.add_option("--ignore-tables",
                      dest    = "ignore_tables",
                      default = None,
                      help    = "The comma-separated list of MySQL tables that should be ignored during migration.")

    parser.add_option("--chunk-size",
                      dest    = "chunksize",
                      default = 1000000,
                      help    = "Maximum number of rows to migrate per transaction (with OFFSET/LIMIT).")


    # GET COMMAND LINE OPTIONS
    (options, arguments) = parser.parse_args()

    if options.DEBUG: logger.setLevel(logging.DEBUG)

    if not options.mysqldb:
        parser.print_help()
        sys.exit(0)

    if not options.mysqlpasswd:
        options.mysqlpasswd = getpass("Input password for MySQL account: ").strip()
    if not options.pgpasswd:
        options.pgpasswd = getpass("Input password for PostgreSQL account: ").strip()

    return parser, options

def get_pg_data_type(mydatatype):
    '''
    '''
    if isinstance(mydatatype, mysql.ENUM):

        # CREATE A STRING COLUMN WITH THE APPROPRIATE LENGTH
        pgcoltype = String(mydatatype.length)

    # MYSQL SET DATA TYPE
    elif isinstance(mydatatype, mysql.SET):

        # CREATE STRING DATA TYPE WITH MAXIMUM SET LENGTH
        # ADD NUMBER OF ITEMS -1 TO ACCOUNT FOR COMMAS
        length = sum(len(v) for v in mydatatype.values) + len(mydatatype.values) - 1
        pgcoltype = String(length)

    # VARIOUS TEXT DATA TYPES
    elif isinstance(mydatatype, mysql.TEXT) or isinstance(mydatatype, mysql.TINYTEXT) \
    or isinstance(mydatatype, mysql.MEDIUMTEXT) or isinstance(mydatatype, mysql.LONGTEXT):

        # CREATE TEXT COLUMN
        pgcoltype = Text()

    # VARIOUS BINARY DATA TYPES
    elif isinstance(mydatatype, mysql.BLOB) or isinstance(mydatatype, mysql.TINYBLOB) \
    or isinstance(mydatatype, mysql.MEDIUMBLOB) or isinstance(mydatatype, mysql.LONGBLOB):

        # CREATE LARGE BINARY COLUMN
        pgcoltype = LargeBinary()

    # POSTGRESQL DOES NOT HAVE YEAR DATA TYPE
    elif isinstance(mydatatype, mysql.YEAR):

        #
        pgcoltype = Date()

    # RELY ON THE SQLALCHEMY TYPE AFFINITY
    else:

        # CREATE CORRESPONDING POSTGRESQL DATA TYPE
        pgcoltype = mydatatype.adapt(mydatatype._type_affinity)

    return pgcoltype

def main():
    '''
    '''
    # GET THE COMMAND LINE OPTIONS
    parser, options = parse_options()

    URL = '{driver}://{user}:{passwd}@{host}:{port}/{db}'

    if not options.mysqlpasswd:
        options.mysqlpasswd = getpass("MySQL password for user %s:" % options.mysqluser)
    mysqlconf = {'driver': 'mysql',
                 'user': options.mysqluser,
                 'passwd': options.mysqlpasswd,
                 'host': options.mysqlhost,
                 'port': options.mysqlport,
                 'db': options.mysqldb
                 }

    if not options.pgpasswd:
        options.pgpasswd = getpass("PostgreSQL password for user %s:" % options.pguser)
    pgsqlconf = {'driver': 'postgresql+psycopg2',
                 'user': options.pguser,
                 'passwd': options.pgpasswd,
                 'host': options.pghost,
                 'port': options.pgport,
                 'db': options.pgdb
                 }

    mysqlengine = create_engine(URL.format(**mysqlconf))
    pgsqlengine = create_engine(URL.format(**pgsqlconf))

    # GET THE DATABASE METADATA FROM MYSQL
    mymetadata = MetaData(bind=mysqlengine)
    mymetadata.reflect()

    # NEW POSTGRESQL METADATA
    pgmetadata = MetaData(bind=pgsqlengine)

    # REMOVE TABLES THAT SHOULD NOT BE MIGRATED FROM METADATA
    if options.ignore_tables:
        for tablename in options.ignore_tables.split(','):
            table = mymetadata.tables[tablename]
            mymetadata.remove(table)

    # ONLY INCLUDE SPECIFIED TABLES
    elif options.tables:
        tables_to_keep = options.tables.split(',')

        for tablename in mymetadata.tables.keys():
            if tablename not in tables_to_keep:
                table = mymetadata.tables[tablename]
                mymetadata.remove(table)

    # KEEP TRACK OF TABLE NAMES BETWEEN METADATA WHICH IS REQUIRED FOR DATA MIGRATION
    tablemapping = {}
    datetime_tables = []

    # MIGRATE ALL MYSQL TABLES
    for tablename in mymetadata.tables:

        # GET MYSQL TABLE
        mytable = mymetadata.tables[tablename]

        # CREATE NEW POSTGRESQL TABLE WITH LOWER CASE NAME IN SPECIFIED SCHEMA
        pgtable = Table(convert(mytable.name), pgmetadata, schema=options.pgschema)

        # KEEP TRACK OF TABLE MAPPING BETWEEN DATABASES
        tablemapping[mytable.name] = '{table.schema}.{table.name}'.format(table=pgtable)

        # MIGRATE MYSQL COLUMNS
        for mycolumn in mytable.columns:

            # GET MYSQL COLUMN DATA TYPE
            mydatatype = mycolumn.type

            # USE THIS TO GET CORRESPONDING POSTGRESQL DATA TYPE
            pgdatatype = get_pg_data_type(mydatatype)

            # OVERRIDE NULLABLE FLAG TO ALLOW POSTGRESQL NULL VALUES INSTEAD OF 0000-00-00 00:00:00
            if isinstance(mydatatype, mysql.DATETIME):
                mycolumn.nullable = True
                datetime_tables.append(tablename)

            # CREATE POSTGRESQL COLUMN
            pgcolumn = Column(mycolumn.name.lower(), pgdatatype,
                              autoincrement=mycolumn.autoincrement,
                              nullable=mycolumn.nullable,
                              primary_key=mycolumn.primary_key,
                              unique=mycolumn.unique)

            # ADD POSTGRESQL COLUMN TO TABLE
            pgtable.append_column(pgcolumn)

        # MIGRATE TABLE INDEXES
        for myindex in mytable.indexes:

            # GET THE POSTGRESQL COLUMNS OF THE INDEX
            pgcolumns = [pgtable.c[column.name.lower()] for column in myindex.columns]

            # CREATE A UNIQUE INDEX NAME / TWO COLUMNS SHOULD BE ENOUGH
            name = '_'.join(('idx', pgtable.name, '_'.join(c.name for c in pgcolumns[:2])))

            # POSTGRESQL IDENTIFIERS MUST NOT EXCEED 63 CHARACTERS
            if len(name) > 63:

                # USE FIRST COLUMN AS INDEX NAME
                name = '_'.join(('idx', pgtable.name, '_'.join(c.name for c in pgcolumns[:1])))

            # CREATE INDEX IN POSTGRESQL TABLE
            pgindex = Index(name, *pgcolumns, unique=myindex.unique, table=pgtable)

        # DEBUG TABLE MIGRATION INFO
        logger.debug("Migrated MySQL table schema {0} to PostgreSQL {1}.{2}".format(mytable.name, pgtable.schema, pgtable.name))

    # CREATE THE DATABASE IN POSTGRESQL
    pgmetadata.drop_all(checkfirst=True)
    pgmetadata.create_all()

    # COPY EACH TABLE THROUGH PIPE
    if options.data:
        pgenviron = os.environ.copy()
        pgenviron['PGPASSWORD'] = options.pgpasswd

        # ITERATE THROUGH MAPPED TABLES AND PIPE DATA
        for mytable, pgtable in sorted(tablemapping.items()):

            # GET TOTAL NUMBER OF ROWS IN MYSQL TABLE TO CREATE CHUNKS
            rowcount = mysqlengine.execute("SELECT COUNT(*) FROM {0}.{1}".format(options.mysqldb, mytable)).scalar()

            logger.debug('{0:,} rows found in table {1}'.format(rowcount,mytable))

            # SPLIT MYSQL RESULT SET WITH LIMIT AND OFFSET
            for offset in xrange(0, rowcount, options.chunksize):

                logger.debug("\tCopying rows {0:,} to {1:,} from MySQL to PostgreSQL".format(offset, offset+options.chunksize))

                # COMPILE MYSQL COMMAND
                mysqldump = Popen(['mysql', '-u', options.mysqluser, '-p%s' % options.mysqlpasswd, '-h', options.mysqlhost,
                                   '-D', options.mysqldb, '--show-warnings', '-s', '--compress',
                                   '--connect_timeout=15', '-N', '-e',
                                   'SELECT * FROM {0} LIMIT {1} OFFSET {2}'.format(mytable, options.chunksize, offset)],
                                   stdout=PIPE)

                # POSTGRESQL COPY
                pgcopyargs = ['psql', '-h', options.pghost, '-d', options.pgdb,
                              '-p', str(options.pgport), '-c', "COPY %s FROM STDIN WITH NULL AS 'NULL'" % pgtable]

                # CHECK IF MYSQL TABLE HAS DATETIME AND POSSIBLE JESUS BIRTHDAY
                if mytable in datetime_tables:

                    # OPEN ANOTHER PIPE
                    sed = Popen(['sed','s/0000-00-00 00:00:00/NULL/g'], stdin=mysqldump.stdout, stdout=PIPE)

                    # COPY POSTGRESQL TABLE FROM SED STDOUT
                    pgcopy = Popen(pgcopyargs, stdin=sed.stdout, env=pgenviron)

                    # RETRY UNTIL DATA WAS SUCCESSFULLY COPIED
                    retcode = mysqldump.wait()

                    sed.wait()
                    pgcopy.wait()

                # PROCEED AS NORMAL
                else:

                    # COPY POSTGRESQL TABLE FROM MYSQL STDOUT
                    pgcopy = Popen(pgcopyargs, stdin=mysqldump.stdout, env=pgenviron)

                    retcode = mysqldump.wait()

                    pgcopy.wait()

            logger.debug("Migrated data from MySQL {0}.{1} to PostgreSQL {2}.".format(options.mysqldb, mytable, pgtable))
main()
