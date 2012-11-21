"""
Contains the definitions of the chain table as well as the subset tables for each
polymer type.
"""

from sqlalchemy import Boolean, Column, Index, Integer, String, Table, Text, DefaultClause
from sqlalchemy.dialects.postgresql import ARRAY, DOUBLE_PRECISION

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

chains = Table('chains', metadata,
               Column('chain_id', Integer, primary_key=True),
               Column('biomolecule_id', Integer, nullable=False),
               Column('pdb_chain_id', String(1), nullable=False), # CASE-SENSITIVE BY DEFAULT
               Column('pdb_chain_asu_id', String(1), nullable=False),
               Column('path', PTree, nullable=False),
               Column('chain_type', String(50)),
               Column('title', Text),
               Column('chain_length', Integer),
               Column('chain_seq', Text),
               Column('chain_seq_md5', String(32)),
               Column('rotation', ARRAY(DOUBLE_PRECISION)),
               Column('translation', ARRAY(DOUBLE_PRECISION)),
               Column('is_at_identity', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
               Column('has_disordered_regions', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
               schema=schema)

Index('idx_chains_biomolecule_id_pdb_chain_id', chains.c.biomolecule_id, chains.c.pdb_chain_id, unique=True)
Index('idx_chains_biomolecule_id_pdb_chain_asu_id', chains.c.biomolecule_id, chains.c.pdb_chain_asu_id)
Index('idx_chains_path', chains.c.path, postgresql_using='gist')
Index('idx_chains_chain_seq_md5', chains.c.chain_seq_md5)

chain_comments = {
    "table": "Contains ALL chains found in a PDB structure if at least one residue interacts with another entity.",
    "columns":
    {
        "chain_id": "Primary key.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "pdb_chain_id": "PDB chain identifier.",
        "pdb_chain_asu_id": "The PDB chain identifier this chain originated from.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID.",
        "chain_type": "Polymer type of the chain - derived from mmcif.entity_poly.type",
        "title": "A description of the Chain, with the name of the Chain in parenthesis. Maps to PDB compound name (entity.pdbx_description)",
        "chain_length": "Number of residues in the polymer chain.",
        "chain_seq": "Sequence of the polymer.",
        "chain_seq_md5": "MD5 Hash of the chain sequence.",
        "rotation": "Rotation matrix.",
        "translation": "Translation vector.",
        "is_at_identity": "True if the chain is at identity, i.e. no transformation was performed.",
        "has_disordered_regions": "True if the chain contains at least one disordered region (unobserved residues)."
    }
}

comment_on_table_elements(chains, chain_comments)


## chain subsets: polypeptides, oligonucleotides and polysaccharides


# polypeptides are chains where the type in mmCIF is polypeptide
polypeptides = Table('polypeptides', metadata,
                     Column('chain_id', Integer, primary_key=True, nullable=False),
                     Column('is_wildtype', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_human', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_enzyme', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_drug_target', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_in_uniprot', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_kinase', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     schema=schema)

polypeptide_comments = {
    "table": "Polypeptides are chains where the type in mmCIF is polypeptide.",
    "columns":
    {
        "chain_id": "Primary key of the parent chain.",
        "is_wildtype": "True if the polypeptide does not contain any mutations (deviations from the UniProt sequence).",
        "is_human": "True if the polypeptide is the product of a human gene.",
        "is_enzyme": "True if the mapped UniProt entry has an Enzyme Code (EC).",
        "is_drug_target": "True if the protein is a drug target in the ChEMBL database.",
        "is_in_uniprot": "True if the polypeptide is part of the UniProt knowledgebase.",
        "is_kinase": "True if the polypeptide is part of the UniProt human & mouse kinase collection."
    }
}

comment_on_table_elements(polypeptides, polypeptide_comments)

# oligonucleotides are chains where the type in mmCIF is polyribonucleotide,
# polydeoxyribonucleotide or a hybrid of both
oligonucleotides = Table('oligonucleotides', metadata,
                         Column('chain_id', Integer, primary_key=True, nullable=False),
                         Column('nucleic_acid_type', String(3)), # DNA/RNA/HYB
                         schema=schema)

oligonucleotide_comments = {
    "table": "Oligonucleotides are chains where the type in mmCIF is polyribonucleotide, polydeoxyribonucleotide or a hybrid of both.",
    "columns":
    {
        "chain_id": "Primary key of the parent chain.",
        "nucleic_acid_type": "Type of Nucleid acid. Either DNA, RNA or HYB (Hybrid)."
    }
}

comment_on_table_elements(oligonucleotides, oligonucleotide_comments)

#
polysaccharides = Table('polysaccharides', metadata,
                        Column('chain_id', Integer, primary_key=True, nullable=False),
                        schema=schema)

polysaccharide_comments = {
    "table": "",
    "columns":
    {
        "chain_id": "Primary key of the parent chain."
    }
}

comment_on_table_elements(polysaccharides, polysaccharide_comments)