"""
"""

from sqlalchemy import Column, Index, Integer, String, Table

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import comment_on_table_elements

xrefs = Table('xrefs', metadata,
              Column('xref_id', Integer, primary_key=True),
              Column('entity_type', String(12), nullable=False),
              Column('entity_id', Integer, nullable=False),
              Column('source', String(32), nullable=False),
              Column('xref', String(64), nullable=False),
              Column('description', Text),
              schema=schema)

Index('idx_xrefs_entity_type', xrefs.c.entity_type, xrefs.c.entity_id)
Index('idx_xrefs_xref', xrefs.c.source, xrefs.c.xref)

comments = {
    "table": "Contains cross references between entities in CREDO and external databases.",
    "columns":
    {
        "xref_id": "Primary key of the cross reference.",
        "entity_id": "Primary key of the CREDO entity.",
        "entity_type": "Type of the CREDO entity, e.g. Ligand.",
        "source": "Name of the external database.",
        "xref": "cross reference in the external database.",
        "description": "Further description of the cross reference."
    }
}

comment_on_table_elements(xrefs, comments)