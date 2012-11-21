TRUNCATE TABLE credo.xrefs;

-- PUBMED IDS FOR PDB STRUCTURES
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Structure', s.structure_id, 'PubMed',c.pdbx_database_id_PubMed, c.journal_id_ASTM
FROM        credo.structures s
JOIN        mmcif.citation c ON s.pdb = c.structure_id
WHERE       pdbx_database_id_PubMed > 0
ORDER BY    s.structure_id;

-- UNIPROT ENTRIES FOR CHAINS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref)
SELECT      DISTINCT 'Chain', c.chain_id, 'UniProt', m.uniprot
FROM        credo.chains c
JOIN        credo.biomolecules b USING(biomolecule_id)
JOIN        credo.structures s USING(structure_id)
JOIN        pdb.res_map m
            ON m.pdb = s.pdb
            AND m.pdb_chain_id = c.pdb_chain_asu_id
WHERE       m.uniprot IS NOT NULL
ORDER BY    chain_id, uniprot;

-- CROSS REFERENCES FOR PDB CHAINS FROM MSD SIFTS
  INSERT INTO credo.xrefs(entity_type, entity_id, source, xref)
  SELECT DISTINCT 'Chain', c.chain_id, db_source, db_accession_id
    FROM pdb.map_regions mr
    JOIN credo.structures s ON s.pdb = mr.pdb
    JOIN credo.chains c
         ON c.pdb_chain_asu_id = mr.pdb_chain_id
         AND subptree(path,0,1)::text = s.pdb
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
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'ChEMBL Compound', ci.chembl_id as xref, md.pref_name
FROM        chembl.compound_smiles cs
JOIN        pdbchem.chem_comps cp ON cp.ism = cs.ism
JOIN        chembl.chembl_id_lookup ci ON ci.entity_type = 'COMPOUND' and ci.entity_id = cs.molregno
LEFT JOIN   chembl.molecule_dictionary md on md.molregno = cs.molregno
ORDER BY    2,4;

-- KEGG COMPOUNDS / REQUIRES SCIFDW FOREIGN DATA WRAPPER!!!
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref)
SELECT      DISTINCT 'ChemComp', cp.chem_comp_id, 'KEGG Compound', u.compound_id as xref
FROM        pdbchem.chem_comps cp
JOIN        credo.unichem_pdb_to_kegg u ON cp.het_id = u.het_id
ORDER BY    2,4;

-- CHEMBL TARGETS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Chain', x.entity_id, 'ChEMBL Target', td.chembl_id as xref, td.pref_name
FROM        credo.xrefs x
JOIN        chembl.target_dictionary td ON td.protein_accession = x.xref
WHERE       x.source = 'UniProt'
ORDER BY    2,4;

-- CHEMBL PROTEIN THERAPEUTICS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Chain', x.entity_id, 'ChEMBL Protein Therapeutic', td.chembl_id as xref, td.pref_name
FROM        credo.xrefs x
JOIN        chembl.target_dictionary td ON td.protein_accession = x.xref
JOIN        chembl.protein_therapeutics p ON p.protein_sequence = td.protein_sequence
WHERE       x.source = 'UniProt'
ORDER BY    2,4;

-- CHEMBL DOCS FOR PDB STRUCTURES
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Structure', x.entity_id, 'ChEMBL Document', ci.chembl_id as xref, d.title
FROM        credo.xrefs x
JOIN        chembl.docs d ON x.xref::int = d.pubmed_id
JOIN        chembl.chembl_id_lookup ci ON ci.entity_type = 'DOCUMENT' and ci.entity_id = d.doc_id
WHERE       x.source = 'PubMed'
ORDER BY    2,4;

-- CHEMBL ASSAYS FOR STRUCTURES (SAME PUBLICATION) AND CHAINS (SAME PROTEIN ACCESSION)
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Chain', x2.entity_id, 'ChEMBL Assay', a.chembl_id as xref, a.description
FROM        credo.xrefs x1
            -- SAME PULICATION
JOIN        chembl.docs d ON x1.source = 'PubMed' AND x1.xref::int = d.pubmed_id
JOIN        credo.biomolecules b ON b.structure_id = x1.entity_id
JOIN        credo.chains c ON c.biomolecule_id = b.biomolecule_id
            -- SAME PROTEIN ACCESSION
JOIN        credo.xrefs x2 ON x2.entity_type = 'Chain' and x2.entity_id = c.chain_id
JOIN        chembl.target_dictionary td ON x2.source = 'UniProt' AND x2.xref = td.protein_accession
            -- GET ASSAYS
JOIN        chembl.assay2target at ON at.tid = td.tid
JOIN        chembl.assays a ON a.assay_id = at.assay_id
ORDER BY    2,4;

-- CHEMBL ACTIVITIES FOR LIGANDS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
WITH        assays AS
            (
            SELECT  DISTINCT x.entity_id as chain_id, a.assay_id, a.description
            FROM    credo.xrefs x
            JOIN    chembl.assays a ON x.source = 'ChEMBL Assay' AND x.xref = a.chembl_id
            ),
            ligands AS
            (
            SELECT  DISTINCT l.ligand_id, i.entity_id AS molregno
            FROM    credo.xrefs x
            JOIN    chembl.chembl_id_lookup i ON x.source = 'ChEMBL Compound' AND x.xref = i.chembl_id
            JOIN    pdbchem.chem_comps cp ON cp.chem_comp_id = x.entity_id
            JOIN    credo.ligands l ON l.ligand_name = cp.het_id
            )
SELECT      DISTINCT 'Ligand', l.ligand_id, 'ChEMBL Activity', a.activity_id,
            a.standard_type || ' ' || a.relation || ' ' || a.standard_value || ' ' || a.standard_units || COALESCE(' (' || y.description || ')','')
FROM        chembl.activities a
JOIN        assays y ON y.assay_id = a.assay_id
JOIN        ligands l ON l.molregno = a.molregno
            -- LINK LIGANDS AND CHAINS
JOIN        credo.binding_site_residues b ON b.ligand_id = l.ligand_id
JOIN        credo.residues r ON r.residue_id = b.residue_id AND r.chain_id = y.chain_id
ORDER BY    2,4;

-- CHEMBL INFERRED ACTIVITIES FOR LIGANDS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Ligand', l.ligand_id, 'ChEMBL Inferred Activity', a.activity_id,
            a.standard_type || ' ' || a.relation || ' ' || a.standard_value || ' ' || a.standard_units || COALESCE(' (' || ay.description || ')','')
FROM        chembl.activities a
JOIN        chembl.assays ay USING(assay_id)
JOIN        chembl.assay2target att USING(assay_id)
JOIN        chembl.target_dictionary td USING(tid)
            -- TARGETS
JOIN        credo.xrefs xr1 ON xr1.source = 'UniProt' AND xr1.xref = td.protein_accession
            -- COMPOUNDS
JOIN        chembl.compound_smiles cs ON cs.molregno = a.molregno
JOIN        pdbchem.chem_comps cp ON cp.ism = cs.ism
            -- LINK TO CREDO
JOIN        credo.chains c ON c.chain_id = xr1.entity_id AND xr1.entity_type = 'Chain'
JOIN        credo.residues r ON r.chain_id = c.chain_id
JOIN        credo.binding_site_residues b ON b.residue_id = r.residue_id
JOIN        credo.ligands l ON l.ligand_id = b.ligand_id AND l.ligand_name = cp.het_id
WHERE       ay.assay_type = 'B'
            AND standard_units = 'nM'
            AND att.confidence_score >= 8
ORDER BY    2;

-- ENSEMBL VARIATIONS
INSERT      INTO credo.xrefs(entity_type, entity_id, source, xref, description)
SELECT      DISTINCT 'Residue', p.residue_id, 'EnsEMBL Variation', v.variation_id, v.variation_name
FROM        credo.peptides p
JOIN        variations.variation_to_pdb vp USING(res_map_id)
JOIN        variations.variation_to_uniprot vu using(variation_to_uniprot_id)
JOIN        variations.variations v USING(variation_id)
ORDER BY    2,4;

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
