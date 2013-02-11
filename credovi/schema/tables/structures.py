from sqlalchemy import Column, Date, Float, Index, Integer, String, Table, Text
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.dialects.postgresql import ARRAY, VARCHAR

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import comment_on_table_elements

structures = Table('structures', metadata,
                    Column('structure_id', Integer, nullable=False),
                    Column('pdb', String(4), nullable=False),
                    Column('title', Text),
                    Column('authors', Text),
                    Column('method', String(32)),
                    Column('deposition_date', Date),
                    Column('modified_date', Date),
                    Column('resolution', Float(3,1)),
                    Column('r_factor', Float(4,3)),
                    Column('r_free', Float(3,3)),
                    Column('ph', Float(2,1)),
                    Column('dpi_theoretical_min', Float(3,2)),
                    Column('dpi', Float(3,2)),
                    Column('num_biomolecules', Integer),
                    schema=schema)

PrimaryKeyConstraint(structures.c.structure_id, deferrable=True, initially='deferred')
Index('idx_structures_pdb', structures.c.pdb, unique=True)

comments = {
    "table": "Represents a structure from a PDB entry.",
    "columns":
    {
        "structure_id": "Primary key.",
        "pdb": "PDB 4-letter code.",
        "title": "Title of the PDB deposition.",
        "authors": "Concatenated string of author names.",
        "method": "Method used in the experiment (exptl.method).",
        "deposition_date": "Date the entry first entered the PDB database in the form yyyy-mm-dd. Taken from the PDB HEADER record (database_PDB_rev.date_original).",
        "modified_date": "Date the latest PDB revision took place. Taken from the REVDAT record (database_PDB_rev.date).",
        "resolution": "Smallest value for the interplanar spacings for the reflection dataused in the refinement in angstroms. This is called the highest resolution (refine.ls_d_res_high).",
        "r_factor": "Residual factor R for reflections* that satisfy the reflns.observed_criterion were included in the refinement (when the refinement included the calculation of a free R factor) (refine.ls_R_factor_R_work).",
        "r_free": "Residual factor R for reflections that satisfy the reflns.observed_criterion that were used as test reflections (i.e. were excluded from the refinement) (refine.ls_R_factor_R_free).",
        "ph": "pH at which the crystal was grown (exptl_crystal_grow.pH).",
        "dpi": "Diffraction-component Precision Index.",
        "dpi_theoretical_min": "Theoretical minimum of the DPI.",
        "num_biomolecules": "Number of biological assemblies in the asymmetric unit.",
    }
}

comment_on_table_elements(structures, comments)
