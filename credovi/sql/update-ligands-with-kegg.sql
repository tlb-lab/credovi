/* This script updates ligands with substrate and product information from KEGG */

-- CREATE TEMPORARY TABLES TO SPEED UP THE UPDATE
CREATE TEMP TABLE ligand_to_compound_id AS
SELECT DISTINCT ligand_id, xref as compound_id
  FROM credo.xrefs xl
  JOIN pdbchem.chem_comps cc
       ON xl.entity_type= 'ChemComp' AND xl.entity_id = cc.chem_comp_id
  JOIN credo.ligands l
       ON l.ligand_name = cc.het_id
 WHERE xl.source = 'KEGG Compound';

CREATE INDEX idx_ligand_to_compound_id_compound_id
  ON ligand_to_compound_id USING btree (compound_id) WITH (FILLFACTOR=100);

CREATE INDEX idx_ligand_to_compound_id_ligand_id
  ON ligand_to_compound_id USING btree (ligand_id) WITH (FILLFACTOR=100);


-- MAPPING BETWEEN CHAIN IDS AND EC
CREATE TEMP TABLE chain_id_to_ec AS
SELECT DISTINCT entity_id as chain_id, xref as ec
  FROM credo.xrefs xc
 WHERE xc.source = 'EC' ;

CREATE INDEX idx_chain_id_to_ec_ec
  ON chain_id_to_ec USING btree (ec) WITH (FILLFACTOR=100);

CREATE INDEX idx_chain_id_to_ec_chain_id
  ON chain_id_to_ec USING btree (chain_id) WITH (FILLFACTOR=100);


-- MAPPING BETWEEN LIGANDS AND THE CHAINS THEY INTERACT WITH
CREATE TEMP TABLE ligand_to_prox_chain_id AS
SELECT DISTINCT ligand_id, p.chain_id
  FROM credo.binding_site_residues bsr
  JOIN credo.peptides p ON p.residue_id = bsr.residue_id;

CREATE INDEX ligand_to_prox_chain_id_ligand_id
  ON ligand_to_prox_chain_id USING btree (ligand_id) WITH (FILLFACTOR=100);

CREATE INDEX ligand_to_prox_chain_id_chain_id
  ON ligand_to_prox_chain_id USING btree (chain_id) WITH (FILLFACTOR=100);

-- UPDATE LIGAND SUBSTRATES
UPDATE credo.ligands l
   SET is_substrate = true
  FROM (
        SELECT DISTINCT lc.ligand_id
          FROM kegg.enzyme_substrates es
          JOIN ligand_to_compound_id lc ON es.compound_id = lc.compound_id
          JOIN chain_id_to_ec ce ON ce.ec = es.ec
          JOIN ligand_to_prox_chain_id lpc ON lpc.ligand_id = lc.ligand_id AND lpc.chain_id = ce.chain_id
       ) sq
 WHERE sq.ligand_id = l.ligand_id;

-- UPDATE LIGAND PRODUCTS
UPDATE credo.ligands l
   SET is_product = true
  FROM (
        SELECT DISTINCT lc.ligand_id
          FROM kegg.enzyme_products ep
          JOIN ligand_to_compound_id lc ON ep.compound_id = lc.compound_id
          JOIN chain_id_to_ec ce ON ce.ec = ep.ec
          JOIN ligand_to_prox_chain_id lpc ON lpc.ligand_id = lc.ligand_id AND lpc.chain_id = ce.chain_id
       ) sq
 WHERE sq.ligand_id = l.ligand_id;