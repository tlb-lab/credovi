"""
"""

from sqlalchemy import Column, DefaultClause, Float, Index, Integer, String, Table, Text
from sqlalchemy.sql import func
from sqlalchemy.schema import PrimaryKeyConstraint

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

prot_fragments = Table('prot_fragments', metadata,
                       Column('prot_fragment_id', Integer, nullable=False),
                       Column('biomolecule_id', Integer, nullable=False),
                       Column('chain_id', Integer, nullable=False),
                       Column('path', PTree),
                       Column('sstruct_serial', Integer),
                       Column('sstruct', String(1)),
                       Column('fragment_size', Integer, nullable=False),
                       Column('fragment_seq', Text),
                       Column('prot_fragment_nterm_id', Integer),
                       Column('prot_fragment_cterm_id', Integer),
                       Column('completeness', Float(3,2), DefaultClause('0'), nullable=False),
                       schema=schema)

PrimaryKeyConstraint(prot_fragments.c.prot_fragment_id, deferrable=True, initially='deferred')
Index('idx_prot_fragments_biomolecule_id', prot_fragments.c.biomolecule_id)
Index('idx_prot_fragments_chain', prot_fragments.c.chain_id, prot_fragments.c.sstruct_serial, unique=True)
Index('idx_prot_fragments_path', prot_fragments.c.path, postgresql_using='gist')
Index('idx_prot_fragments_seq', prot_fragments.c.fragment_seq)

prot_fragment_comments = {
    "table": "Contains information about secondary structure fragments in proteins.",
    "columns":
    {
        "prot_fragment_id": "Primary key of the protein fragment.",
        "biomolecule_id": "Primary key of the biomolecule the secondary structure fragment belongs to.",
        "chain_id": "Primary key of the chain the secondary structure belongs to.",
        "path": "ptree",
        "sstruct_serial": "Secondary structure serial number inside the chain.",
        "sstruct": "Secondary structure DSSP type.",
        "fragment_size": "Number of residues that form the protein fragment.",
        "fragment_seq": "Fragment one-letter-code sequence.",
        "prot_fragment_nterm_id": "N-terminal protein fragment neighbour.",
        "prot_fragment_cterm_id": "C-terminal protein fragment neighbour.",
        "completeness": "Simple measure to quantify the completeness of the protein fragment: is is the number of its residues that are in CREDO divided by the size of the fragment."
    }
}

comment_on_table_elements(prot_fragments, prot_fragment_comments)

prot_fragment_residues = Table('prot_fragment_residues', metadata,
                               Column('prot_fragment_id', Integer, nullable=False),
                               Column('residue_id', Integer, nullable=False),
                               schema=schema)

PrimaryKeyConstraint(prot_fragment_residues.c.prot_fragment_id,
                     prot_fragment_residues.c.residue_id, deferrable=True,
                     initially='deferred')
Index('idx_prot_fragment_residues_residue_id', prot_fragment_residues.c.residue_id)

prot_fragment_residue_comments = {
    "table": "Mapping between protein fragments and residues in CREDO - The database contains only residue interacting with other entities though, so some residues are missing.",
    "columns":
    {
        "prot_fragment_id": "Primary key of the protein fragment.",
        "residue_id": "Primary key of the residue that is part of the protein fragment."
    }
}

comment_on_table_elements(prot_fragment_residues, prot_fragment_residue_comments)