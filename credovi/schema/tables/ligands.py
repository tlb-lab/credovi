from sqlalchemy import (Boolean, Column, DDL, DefaultClause, Float, Index, Integer,
                        LargeBinary, String, Table, Text)
from sqlalchemy.event import listen
from sqlalchemy.schema import PrimaryKeyConstraint
from sqlalchemy.dialects.postgresql import ARRAY, REAL

from credovi.schema import metadata, schema
from credovi.util.sqlalchemy import ArrayXi, Cube, PTree, comment_on_table_elements

ligands = Table('ligands', metadata,
                Column('ligand_id', Integer),
                Column('biomolecule_id', Integer, nullable=False),
                Column('path', PTree, nullable=False),
                Column('entity_serial', Integer, nullable=False),
                Column('pdb_chain_id', String(1), nullable=False),
                Column('ligand_name', String(64), nullable=False),
                Column('res_num', Integer),
                Column('num_hvy_atoms', Integer),
                Column('ism', Text),
                Column('gini_index_contacts', Float(4,3)),
                Column('is_at_identity', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_incomplete', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_disordered', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_clashing', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_substrate', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_product', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_cofactor', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_promiscuous', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                Column('is_drug_target_int', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                schema=schema)

PrimaryKeyConstraint(ligands.c.ligand_id, deferrable=True, initially='deferred')
Index('idx_ligands_biomolecule_id', ligands.c.biomolecule_id, ligands.c.entity_serial, unique=True)
Index('idx_ligands_pdb', ligands.c.biomolecule_id, ligands.c.pdb_chain_id, ligands.c.res_num, ligands.c.ligand_name, unique=True)
Index('idx_ligands_name', ligands.c.ligand_name, ligands.c.pdb_chain_id, ligands.c.res_num)
Index('idx_ligands_path', ligands.c.path, postgresql_using='gist')

ligand_comments = {
    "table": "A ligand is a set of up to 10 residues.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "path": "ptree path in the form PDB/assembly serial number/PDB chain ID[/PDB residue name`PDB residue number].",
        "entity_serial": "Serial number of the ligand as entity. Only used internally for inserts, updates, etc.",
        "pdb_chain_id": "PDB Chain identifier of the ligand.",
        "ligand_name": "HET-ID or concatenation of them in case of polymer ligands.",
        "res_num": "PDB residue number of the ligand if single residue, otherwise NULL.",
        "num_hvy_atoms": "Number of heavy atoms of the ligand IN the PDB entry.",
        "ism": "Isomeric SMILES.",
        "gini_index_contacts": "Gini index for the contacts that the ligand atom form with other atoms (PRO,RNA,DNA).",
        "is_at_identity": "True if the ligand is at identity, i.e. no transformation was performed.",
        "is_incomplete": "True if the ligand has less atoms than the idealised version from PDBeChem.",
        "is_disordered": "True if at least one ligand atom is disordered.",
        "is_clashing": "True if the ligand is clashing with other residues.",
        "is_substrate": "",
        "is_product": "",
        "is_cofactor": "",
        "is_promiscuous": "",
        "is_drug_target_int": ""
    }
}

comment_on_table_elements(ligands, ligand_comments)

ligand_components = Table('ligand_components', metadata,
                          Column('ligand_component_id', Integer),
                          Column('ligand_id', Integer, nullable=False),
                          Column('residue_id', Integer, nullable=False),
                          Column('het_id', String(3), nullable=False),
                          schema=schema)

PrimaryKeyConstraint(ligand_components.c.ligand_component_id, deferrable=True, initially='deferred')
Index('idx_ligand_components_ligand_id', ligand_components.c.ligand_id, ligand_components.c.residue_id, unique=True)
Index('idx_ligand_components_residue_id', ligand_components.c.residue_id, ligand_components.c.ligand_id, unique=True)

ligand_component_comments = {
    "table": "Contains the residues that a given ligand consists of.",
    "columns":
    {
        "ligand_component_id": "Primary key of the ligand component.",
        "ligand_id": "Primary key of the parent ligand.",
        "residue_id": "Primary key of the residue.",
        "het_id": "Chemical component identifier."
    }
}

comment_on_table_elements(ligand_components, ligand_component_comments)

ligand_fragments = Table('ligand_fragments', metadata,
                         Column('ligand_fragment_id', Integer, nullable=False),
                         Column('biomolecule_id', Integer, nullable=False),
                         Column('ligand_id', Integer, nullable=False),
                         Column('ligand_component_id', Integer, nullable=False),
                         Column('fragment_id', Integer, nullable=False),
                         Column('hit', Integer, nullable=False),
                         Column('num_int_atoms', Integer, DefaultClause('0'), nullable=False),
                         Column('num_contacts', Integer, DefaultClause('0'), nullable=False),
                         Column('num_covalent', Integer, DefaultClause('0'), nullable=False),
                         Column('num_vdw_clash', Integer, DefaultClause('0'), nullable=False),
                         Column('num_vdw', Integer, DefaultClause('0'), nullable=False),
                         Column('num_proximal', Integer, DefaultClause('0'), nullable=False),
                         Column('num_hbond', Integer, DefaultClause('0'), nullable=False),
                         Column('num_weak_hbond', Integer, DefaultClause('0'), nullable=False),
                         Column('num_xbond', Integer, DefaultClause('0'), nullable=False),
                         Column('num_ionic', Integer, DefaultClause('0'), nullable=False),
                         Column('num_metal_complex', Integer, DefaultClause('0'), nullable=False),
                         Column('num_aromatic', Integer, DefaultClause('0'), nullable=False),
                         Column('num_hydrophobic', Integer, DefaultClause('0'), nullable=False),
                         Column('num_carbonyl', Integer, DefaultClause('0'), nullable=False),
                         Column('fcd', Float(5,3), DefaultClause('0'), nullable=False),
                         Column('fad', Float(5,3), DefaultClause('0'), nullable=False),
                         Column('npr1', Float(4,3)),
                         Column('npr2', Float(4,3)),
                         Column('asa_buriedness', Float(4,3)),
                         Column('rel_asa_buriedness', Float(4,3)),
                         Column('aasa_buriedness', Float(4,3)),
                         Column('pasa_buriedness', Float(4,3)),
                         Column('is_root', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                         Column('is_interacting', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                         schema=schema)

PrimaryKeyConstraint(ligand_fragments.c.ligand_fragment_id, deferrable=True, initially='deferred')
Index('idx_ligand_fragments_biomolecule_id', ligand_fragments.c.biomolecule_id)
Index('idx_ligand_fragments_ligand_id', ligand_fragments.c.ligand_id)
Index('idx_ligand_fragments_ligand_component_id', ligand_fragments.c.ligand_component_id)
Index('idx_ligand_fragments_fragment_id', ligand_fragments.c.fragment_id, ligand_fragments.c.hit)

ligand_fragment_comments = {
    "table": "ligand fragments are the result of mapping the RECAP fragments of chemical components onto ligands and their atoms.",
    "columns":
    {
        "ligand_fragment_id": "Primary key of the ligand fragment.",
        "biomolecule_id": "Primary key of the parent biomolecule.",
        "ligand_id": "Primary key of the parent ligand.",
        "ligand_component_id": "Primary key of the parent ligand component.",
        "fragment_id": "Primary key of the fragment that was mapped to create the ligand fragment.",
        "hit": "Hit number of the fragment; sometimes a fragment can be found more than once in a ligand.",
        "asa_buriedness":"Ratio of the accessible surface area of the fragment that is buried in the bound state.",
        "aasa_buriedness":"Ratio of the apolar accessible surface area of the fragment that is buried in the bound state.",
        "pasa_buriedness":"Ratio of the polar accessible surface area of the fragment that is buried in the bound state."
    }
}

comment_on_table_elements(ligand_fragments, ligand_fragment_comments)

ligand_fragment_atoms = Table('ligand_fragment_atoms', metadata,
                              Column('ligand_fragment_atom_id', Integer, nullable=False),
                              Column('ligand_id', Integer, nullable=False),
                              Column('ligand_fragment_id', Integer, nullable=False),
                              Column('atom_id', Integer, nullable=False),
                              schema=schema)

PrimaryKeyConstraint(ligand_fragment_atoms.c.ligand_fragment_atom_id, deferrable=True, initially='deferred')
Index('idx_ligand_fragment_atoms_ligand_fragment_id', ligand_fragment_atoms.c.ligand_fragment_id, ligand_fragment_atoms.c.atom_id, unique=True)
Index('idx_ligand_fragment_atoms_ligand_id', ligand_fragment_atoms.c.ligand_id)
Index('idx_ligand_fragment_atoms_atom_id', ligand_fragment_atoms.c.atom_id)

ligand_molstrings = Table('ligand_molstrings', metadata,
                          Column('ligand_id', Integer, nullable=False, autoincrement=False),
                          Column('ism', Text, nullable=False),
                          Column('pdb', Text, nullable=False),
                          Column('sdf', Text, nullable=False),
                          Column('oeb', LargeBinary, nullable=False),
                          schema=schema)

PrimaryKeyConstraint(ligand_molstrings.c.ligand_id, deferrable=True, initially='deferred')

ligand_molstring_comments = {
    "table": "Contains the chemical structures of all ligands with at least 7 heavy atoms in various formats.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand.",
        "ism": "Isomeric SMILES.",
        "pdb": "PDB format.",
        "sdf": "SDF format.",
        "oeb": "OpenEye binary format."
    }
}

comment_on_table_elements(ligand_molstrings, ligand_molstring_comments)


ligand_usr = Table('ligand_usr', metadata,
                   Column('ligand_id', Integer, nullable=False, autoincrement=False),
                   Column('npr1', Float(4,3)),
                   Column('npr2', Float(4,3)),
                   Column('usr_space', Cube, nullable=False),
                   Column('usr_moments', ARRAY(Float, dimensions=1), nullable=False),
                   schema=schema)

PrimaryKeyConstraint(ligand_usr.c.ligand_id, deferrable=True, initially='deferred')

# create gist index on n-dimensional space
Index ('idx_ligand_usr_usr_space', ligand_usr.c.usr_space, postgresql_using='gist')

listen(ligand_usr, "after_create",
       DDL("ALTER TABLE %(fullname)s CLUSTER ON idx_ligand_usr_usr_space;",
           on='postgresql'))

ligand_usr_comments = {
    "table": "Contains the USR moments of all ligands with at least 7 heavy atoms. The usr_space columns contains the first 12 default moments as CUBE to make use of the GIST index. usr_moments contains the USRCAT moments.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand.",
        "usr_space": "Default USR moments as PostgreSQL CUBE data type to make use of the GIST index and other n-dimensional space functions.",
        "usr_moments": "USRCAT moments."
    }
}

comment_on_table_elements(ligand_usr, ligand_usr_comments)


ligand_eff = Table('ligand_eff', metadata,
                   Column('ligand_id', Integer, primary_key=True, nullable=False, autoincrement=False),
                   Column('activity_id', Integer, primary_key=True, nullable=False, autoincrement=False),
                   Column('assay_chembl_id', Text),
                   Column('assay_description', Text),
                   Column('standard_type', Text),
                   Column('relation', String(50)),
                   Column('standard_value', REAL),
                   Column('standard_units', Text),
                   Column('p', REAL),
                   Column('bei', REAL),
                   Column('sei', REAL),
                   Column('activity_comment', Text),
                   schema=schema)

PrimaryKeyConstraint(ligand_eff.c.ligand_id, ligand_eff.c.activity_id,
                     deferrable=True, initially='deferred')

ligand_eff_comments = {
    "table": "Contains ligand efficiency data from ChEMBL.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand."
    }
}

comment_on_table_elements(ligand_eff, ligand_eff_comments)


hetatms = Table('hetatms', metadata,
                Column('hetatm_id', Integer, nullable=False),
                Column('ligand_id', Integer, nullable=False),
                Column('ligand_component_id', Integer, nullable=False),
                Column('atom_id', Integer, nullable=False),
                schema=schema)

PrimaryKeyConstraint(hetatms.c.hetatm_id, deferrable=True, initially='deferred')
Index('idx_hetatms_ligand_id', hetatms.c.ligand_id)
Index('idx_hetatms_ligand_component_id', hetatms.c.ligand_component_id)
Index('idx_hetatms_atom_id', hetatms.c.atom_id)

hetatm_comments = {
    "table": "Contains the atoms a given ligand consists of.",
    "columns":
    {

        "hetatm_id": "Primary key of the hetatm.",
        "ligand_id": "Primary key of the parent ligand.",
        "ligand_component_id": "Primary key of the parent ligand component.",
        "atom_id": "Primary key of the atom."
    }
}

comment_on_table_elements(hetatms, hetatm_comments)

#
binding_sites = Table('binding_sites', metadata,
                      Column('ligand_id', Integer, autoincrement=False, nullable=False),
                      Column('has_incomplete_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      Column('has_non_std_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      Column('has_mod_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      Column('has_mut_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      Column('has_mapped_var', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      Column('is_kinase', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                      schema=schema)

PrimaryKeyConstraint(binding_sites.c.ligand_id, deferrable=True,
                     initially='deferred')

binding_sites_comments = {
    "table": "Contains more info about the ligand binding sites (protein only!).",
    "columns":
    {
        "ligand_id": "Primary key of the ligand.",
        "cath_dmns": "Array of cath domains that are mapped onto this binding site.",
        "scop_pxs": "Array of SCOP px identifiers that are mapped onto this binding site.",
        "is_mutated": "True if at least one amino acid differs from the canonical UniProt sequence.",
        "has_incomplete_res": "True if at least one amino acid that is part if this binding site has missing atoms.",
        "has_mapped_var": "True if at least one residue can be linked to a variation.",
        "is_kinase": "True if the binding site polypeptide is part of the UniProt human & mouse kinase collection."
    }
}

# link between ligands and the residues they interact with
binding_site_residues = Table('binding_site_residues', metadata,
                              Column('ligand_id', Integer, autoincrement=False, nullable=False),
                              Column('residue_id', Integer,  autoincrement=False, nullable=False),
                              Column('entity_type_bm', Integer, nullable=False),
                              schema=schema)

PrimaryKeyConstraint(binding_site_residues.c.ligand_id, binding_site_residues.c.residue_id,
                     deferrable=True, initially='deferred')
Index('idx_binding_site_residues_residue_id', binding_site_residues.c.residue_id)

binding_site_residues_comments = {
    "table": "Mapping between ligands and the residues they interact with.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand.",
        "residue_id": "Primary key of the residue that is interacting."
    }
}

comment_on_table_elements(binding_site_residues, binding_site_residues_comments)

binding_site_domains = Table('binding_site_domains', metadata,
                              Column('ligand_id', Integer, autoincrement=False, nullable=False),
                              Column('domain_id', Integer, autoincrement=False, nullable=False),
                              schema=schema)

PrimaryKeyConstraint(binding_site_domains.c.ligand_id, binding_site_domains.c.domain_id,
                     deferrable=True, initially='deferred')
Index('idx_binding_site_domains_domain_id', binding_site_domains.c.domain_id,
      binding_site_domains.c.ligand_id, unique=True)

# interactions between ligands
lig_lig_ints = Table('lig_lig_interactions', metadata,
                     Column('lig_lig_interaction_id', Integer, nullable=False, autoincrement=True),
                     Column('biomolecule_id', Integer, nullable=False),
                     Column('lig_bgn_id', Integer, nullable=False),
                     Column('lig_end_id', Integer, nullable=False),
                     Column('path', PTree),
                     Column('is_quaternary', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('is_homo_dimer', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_incomplete_res', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_clash', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_product', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_substrate', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_drug_target_int', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     Column('has_drug_like_ligands', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     schema=schema)

PrimaryKeyConstraint(lig_lig_ints.c.lig_lig_interaction_id, deferrable=True, initially='deferred')
Index('idx_lig_lig_interactions_biomolecule_id', lig_lig_ints.c.biomolecule_id)
Index('idx_lig_lig_interactions_path', lig_lig_ints.c.path, postgresql_using='gist')
Index('idx_lig_lig_interactions_lig_bgn_id', lig_lig_ints.c.lig_bgn_id)
Index('idx_lig_lig_interactions_lig_end_id', lig_lig_ints.c.lig_end_id)

# interactions between ligands and nucleic acids
lig_nuc_ints = Table('lig_nuc_interactions', metadata,
                     Column('lig_nuc_interaction_id', Integer, nullable=False, autoincrement=True),
                     Column('biomolecule_id', Integer, nullable=False),
                     Column('ligand_id', Integer, nullable=False),
                     Column('chain_nuc_id', Integer, nullable=False),
                     Column('path', PTree),
                     Column('is_quaternary', Boolean(create_constraint=False), DefaultClause('false'), nullable=False),
                     schema=schema)

PrimaryKeyConstraint(lig_nuc_ints.c.lig_nuc_interaction_id, deferrable=True, initially='deferred')
Index('idx_lig_nuc_interactions_biomolecule_id', lig_nuc_ints.c.biomolecule_id)
Index('idx_lig_nuc_interactions_path', lig_nuc_ints.c.path, postgresql_using='gist')
Index('idx_lig_nuc_interactions_ligand_id', lig_nuc_ints.c.ligand_id)
Index('idx_lig_nuc_interactions_chain_nuc_id', lig_nuc_ints.c.chain_nuc_id)


binding_site_atom_surface_areas = Table('binding_site_atom_surface_areas', metadata,
                                        Column('ligand_id', Integer, primary_key=True, autoincrement=False),
                                        Column('atom_id', Integer, primary_key=True, autoincrement=False),
                                        Column('asa_apo', Float(4,2)),
                                        Column('asa_bound', Float(4,2)),
                                        Column('asa_delta', Float(4,2)),
                                        schema=schema)

Index('idx_binding_site_atom_surface_areas_atom_id', binding_site_atom_surface_areas.c.atom_id)

binding_site_atom_surface_area_commnents = {
    "table": "Solvent-accessible surface area changes for each atom in a binding site defined by a ligand.",
    "columns":
    {
        "ligand_id": "Primary key of the ligand that defines this binding site.",
        "atom_id": "Primary key of the atom that has a different solvent-exposed surface upon binding - can be a ligand or polymer atom.",
        "asa_apo": "Solvent-accessible surface area of the atom in the apo state.",
        "asa_bound": "Solvent-accessible surface area of the atom in the bound state.",
        "asa_delta": "Change in solvent-accessible surface area between the apo and bound state."
    }
}

comment_on_table_elements(binding_site_atom_surface_areas, binding_site_atom_surface_area_commnents)


# table containing the fuzcav fingerprint for all ligand-defined binding sites
binding_site_fuzcav = Table('binding_site_fuzcav', metadata,
                            Column('ligand_id', Integer, nullable=False, autoincrement=False),
                            Column('calphafp', ArrayXi, nullable=False),
                            Column('repfp', ArrayXi, nullable=False),
                            schema=schema)

PrimaryKeyConstraint(binding_site_fuzcav.c.ligand_id, deferrable=True,
                     initially='deferred')
