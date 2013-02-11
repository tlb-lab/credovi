"""
The grooves table holds protein-nucleic acid complexes.
"""

from sqlalchemy import Column, Index, Integer, String, Table, Boolean, DefaultClause

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import PTree, comment_on_table_elements

grooves = Table('grooves', metadata,
                Column('groove_id', Integer, primary_key=True),
                Column('biomolecule_id', Integer, nullable=False),
                Column('chain_prot_id', Integer, nullable=False),
                Column('chain_nuc_id', Integer, nullable=False),
                Column('path', PTree),
                Column('num_res_prot', Integer, nullable=False),
                Column('num_res_nuc', Integer, nullable=False),
                Column('nucleic_acid_type', String(3), nullable=False), # DNA, RNA, HYB
                Column('has_incomplete_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_quaternary', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                schema=schema)

Index('idx_grooves_chain_prot_id', grooves.c.chain_prot_id)
Index('idx_grooves_chain_nuc_id', grooves.c.chain_nuc_id)
Index('idx_grooves_path', grooves.c.path, postgresql_using='gist')

groove_comments = {
    "table": "Contains all the binary interactions between proteins and nucleotides, henceforth known as grooves.",
    "columns":
    {
        "groove_id": "",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "chain_prot_id": "Primary key of the polypeptide chain.",
        "chain_nuc_id": "Primary key of the oligonucleotide chain.",
        "path": "ptree path in the form PDB/assembly serial number/G:PDB chain prot ID-PDB chain nuc ID",
        "num_res_prot": "Number of polypetide residues that are interacting.",
        "num_res_nuc": "Number of oligonucleotide residues that are interacting.",
        "nucleic_acid_type": "Type of the nucleic acid, either DNA/RNA or hybrid.",
        "has_missing_atoms": "True if at least one residue has missing atoms."
    }
}

comment_on_table_elements(grooves, groove_comments)

groove_residues = Table('groove_residues', metadata,
                        Column('groove_id', Integer, autoincrement=False, primary_key=True),
                        Column('residue_prot_id', Integer, autoincrement=False, primary_key=True),
                        Column('residue_nuc_id', Integer, autoincrement=False, primary_key=True),
                        schema=schema)

Index('idx_groove_residues_residue_bgn_id', groove_residues.c.residue_prot_id)
Index('idx_groove_residues_residue_end_id', groove_residues.c.residue_nuc_id)

groove_residue_comments = {
    "table": "contains the combinations of all residues in an interface that are interacting.",
    "columns":
    {
        "groove_id": "Primary key of the polypeptide-nucleic acid interaction.",
        "residue_prot_id": "Primary key of the polypeptide residue.",
        "residue_nuc_id": "Primary key of the oligonucleotide residue."
    }
}

comment_on_table_elements(groove_residues, groove_residue_comments)
