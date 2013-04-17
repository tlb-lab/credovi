TRUNCATE TABLE credo.xrefs;

-- PUBMED IDS FOR PDB STRUCTURES
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Structure', s.structure_id, 'PubMed',c.pdbx_database_id_PubMed, c.journal_id_ASTM
FROM        credo.structures s
JOIN        mmcif.citation c ON s.pdb = c.structure_id
WHERE       pdbx_database_id_PubMed > 0
ORDER BY    s.structure_id;

-- CROSS REFERENCES FOR PDB CHAINS FROM MSD SIFTS
  INSERT INTO credo.xrefs(entity_type, entity_id, source, xref)
  SELECT DISTINCT 'Chain', c.chain_id, db_source, db_accession_id
    FROM pdb.map_regions mr
    JOIN credo.structures s ON s.pdb = mr.pdb
    JOIN credo.biomolecules b ON s.structure_id = b.structure_id
    JOIN credo.chains c
         ON c.biomolecule_id = b.biomolecule_id
         AND c.pdb_chain_asu_id = mr.pdb_chain_id
ORDER BY 2,3,4;

-- DRUGBANK COMPOUNDS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'DrugBank Compound', d.drugbank_id as xref, d.name
FROM        drugbank.smiles ds
JOIN        drugbank.drugs d USING(drugbank_id)
JOIN        pdbchem.chem_comps cp ON cp.ism = ds.ism
ORDER BY    2,4;

-- DRUGBANK TARGETS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref)
SELECT      DISTINCT 'Chain', x.entity_id, 'DrugBank Target', tx.target_id as xref
FROM        credo.xrefs x
JOIN        drugbank.target_xrefs tx ON tx.xref = x.xref
WHERE       x.source = 'UniProt' AND tx.source = 'UniProtKB'
ORDER BY    2,4;

-- CHEMBL COMPOUNDS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'ChEMBL Compound', u.chembl_id as xref, md.pref_name
FROM        scifdw.unichem_pdb_to_chembl u
JOIN        pdbchem.chem_comps cp ON cp.het_id = u.het_id
LEFT JOIN   chembl.molecule_dictionary md on md.chembl_id = u.chembl_id
ORDER BY    2,4;

-- KEGG COMPOUNDS / REQUIRES SCIFDW FOREIGN DATA WRAPPER!!!
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref)
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'KEGG Compound', u.compound_id as xref
FROM        pdbchem.chem_comps cp
JOIN        scifdw.unichem_pdb_to_kegg u ON cp.het_id = u.het_id
ORDER BY    2,4;

-- CHEBI COMPOUNDS / REQUIRES SCIFDW FOREIGN DATA WRAPPER!!!
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref)
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'ChEBI', u.chebi_id as xref
FROM        pdbchem.chem_comps cp
JOIN        scifdw.unichem_pdb_to_chebi u ON cp.het_id = u.het_id
ORDER BY    2,4;

-- Protein targets from ChEMBL - identified via UniProt
  INSERT INTO credo.xrefs(entity_type, entity_id, source, xref, description)
  SELECT DISTINCT 'Chain', x.entity_id, 'ChEMBL Target', td.chembl_id as xref, td.pref_name
    FROM credo.xrefs x
    JOIN chembl.component_sequences cs ON cs.accession = x.xref
    JOIN chembl.target_components tc ON tc.component_id = cs.component_id
    JOIN chembl.target_dictionary td ON td.tid = tc.tid
   WHERE x.source = 'UniProt'
ORDER BY 2,4;

-- insert ChEMBL biotherapeutics using sequence match
  INSERT INTO credo.xrefs(entity_type, entity_id, source, xref, description)
  SELECT DISTINCT 'Chain', c.chain_id, 'ChEMBL Biotherapeutic', md.chembl_id, md.pref_name
    FROM credo.chains c
    JOIN chembl.bio_component_sequences bcs ON c.chain_seq = bcs.sequence
    JOIN chembl.biotherapeutic_components bc ON bc.component_id = bcs.component_id
    JOIN chembl.biotherapeutics b ON b.molregno = bc.molregno
    JOIN chembl.molecule_dictionary md ON md.molregno = b.molregno
ORDER BY 2,4;

-- chembl docs for pdb structures
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Structure', x.entity_id, 'ChEMBL Document', ci.chembl_id as xref, d.title
FROM        credo.xrefs x
JOIN        chembl.docs d ON x.xref::int = d.pubmed_id
JOIN        chembl.chembl_id_lookup ci ON ci.entity_type = 'DOCUMENT' and ci.entity_id = d.doc_id
WHERE       x.source = 'PubMed'
ORDER BY    2,4;

-- chembl assays for structures (same publication) and chains (same protein accession)
  INSERT INTO credo.xrefs(entity_type, entity_id, source, xref, description)
  SELECT DISTINCT 'Chain', x2.entity_id, 'ChEMBL Assay', a.chembl_id as xref, a.description
    FROM credo.xrefs x1
         -- same publication
    JOIN chembl.docs d ON x1.source = 'PubMed' AND x1.xref::int = d.pubmed_id
    JOIN credo.biomolecules b ON b.structure_id = x1.entity_id
    JOIN credo.chains c ON c.biomolecule_id = b.biomolecule_id
         -- same UniProt accession
    JOIN credo.xrefs x2 ON x2.entity_type = 'Chain' and x2.entity_id = c.chain_id
    JOIN chembl.component_sequences cs ON x2.source = 'UniProt' AND x2.xref = cs.accession
         -- get the assays
    JOIN chembl.target_components tc ON tc.component_id = cs.component_id
    JOIN chembl.assays a ON a.tid = tc.tid AND a.doc_id = d.doc_id
ORDER BY 2,4;

-- UPDATE FLAG FOR LIGAND DRUG-TARGET INTERACTIONS
CREATE TEMP TABLE ligand_to_drugbank AS
SELECT DISTINCT l.ligand_id, dd.drugbank_id
  FROM credo.ligands l
  JOIN pdbchem.chem_comps cc ON cc.het_id = l.ligand_name
  JOIN credo.xrefs xr ON xr.entity_id = cc.chem_comp_id
  JOIN drugbank.drugs dd ON dd.drugbank_id = xr.xref
 WHERE xr.entity_type = 'ChemComp'
       AND xr.source = 'DrugBank Compound'
       AND dd.groups @> ARRAY['approved'];

CREATE INDEX idx_ligand_to_drugbank_ligand_id ON ligand_to_drugbank USING btree (ligand_id);
CREATE INDEX idx_ligand_to_drugbank_drugbank_id ON ligand_to_drugbank USING btree (drugbank_id);

UPDATE credo.ligands l
   SET is_drug_target_int = true
  FROM (
        SELECT DISTINCT ld.ligand_id
          FROM ligand_to_drugbank ld
          JOIN credo.binding_site_residues bs ON bs.ligand_id = ld.ligand_id
          JOIN credo.peptides p ON bs.residue_id = p.residue_id
          JOIN credo.xrefs xc ON xc.entity_type = 'Chain' and xc.entity_id = p.chain_id
          JOIN drugbank.drug_to_target dt ON dt.drugbank_id = ld.drugbank_id AND dt.target_id = xc.xref::int
         WHERE xc.source = 'DrugBank Target'
       ) sq
  WHERE l.ligand_id = sq.ligand_id;
