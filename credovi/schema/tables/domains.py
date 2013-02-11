from sqlalchemy import Column, Index, Integer, String, Table, Text
from sqlalchemy.event import listen
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.dialects.postgresql import ARRAY, REAL

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import comment_on_table_elements

domains = Table('domains', metadata,
                Column('domain_id', Integer, autoincrement=False, nullable=False),
                Column('db_source', String(12), nullable=False),
                Column('db_accession_id', String(12), nullable=False),
                Column('description', Text),
                schema=schema)

PrimaryKeyConstraint(domains.c.domain_id, deferrable=True, initially='deferred')
Index('idx_domains_db_accession_id', domains.c.db_source, domains.c.db_accession_id, unique=True)

domain_peptides = Table('domain_peptides', metadata,
                          Column('domain_id', Integer, autoincrement=False, nullable=False),
                          Column('residue_id', Integer, nullable=False),
                          schema=schema)

PrimaryKeyConstraint(domain_peptides.c.domain_id, domain_peptides.c.residue_id,
                     deferrable=True, initially='deferred')
Index('idx_domain_peptides_residue_id', domain_peptides.c.residue_id,
      domain_peptides.c.domain_id, unique=True)
