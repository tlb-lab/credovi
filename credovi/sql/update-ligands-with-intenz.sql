/* This script updates ligands with substrate and product information from IntEnz.
   The IntEnz data adds several hundred more hits to the KEGG data. */

-- mapping between chain ids and ec
CREATE TEMP TABLE chain_id_to_ec AS
SELECT DISTINCT entity_id as chain_id, xref as ec
  FROM credo_dev.xrefs xc
 WHERE xc.source = 'EC' ;

CREATE INDEX idx_chain_id_to_ec_ec
  ON chain_id_to_ec USING btree (ec) WITH (FILLFACTOR=100);

CREATE INDEX idx_chain_id_to_ec_chain_id
  ON chain_id_to_ec USING btree (chain_id) WITH (FILLFACTOR=100);

-- create temporary tables to speed up the update
CREATE TEMP TABLE ligand_to_chebi_id AS
SELECT DISTINCT ligand_id, chebi_id
  FROM scifdw.unichem_pdb_to_chebi u
  JOIN credo_dev.ligands l ON l.ligand_name = u.het_id;

CREATE INDEX idx_ligand_to_chebi_id_chebi_id
  ON ligand_to_chebi_id USING btree (chebi_id) WITH (FILLFACTOR=100);

CREATE INDEX idx_ligand_to_chebi_id_ligand_id
  ON ligand_to_chebi_id USING btree (ligand_id) WITH (FILLFACTOR=100);

-- mapping between ligands and the chains they interact with
CREATE TEMP TABLE ligand_to_prox_chain_id AS
SELECT DISTINCT ligand_id, p.chain_id
  FROM credo_dev.binding_site_residues bsr
  JOIN credo_dev.peptides p ON p.residue_id = bsr.residue_id;

CREATE INDEX ligand_to_prox_chain_id_ligand_id
  ON ligand_to_prox_chain_id USING btree (ligand_id) WITH (FILLFACTOR=100);

CREATE INDEX ligand_to_prox_chain_id_chain_id
  ON ligand_to_prox_chain_id USING btree (chain_id) WITH (FILLFACTOR=100);

-- update ligand substrates
UPDATE credo_dev.ligands l
   SET is_substrate = true
  FROM (
        SELECT DISTINCT lc.ligand_id
          FROM intenz.reactions r, ligand_to_chebi_id lc, chain_id_to_ec ce,
               ligand_to_prox_chain_id lpc
         WHERE r.ec = ce.ec
               AND lpc.ligand_id = lc.ligand_id AND lpc.chain_id = ce.chain_id
               AND lc.chebi_id = ANY(r.reactants)
       ) sq
 WHERE is_substrate = false AND sq.ligand_id = l.ligand_id;

-- update ligand products
UPDATE credo_dev.ligands l
   SET is_product = true
  FROM (
        SELECT DISTINCT lc.ligand_id
          FROM intenz.reactions r, ligand_to_chebi_id lc, chain_id_to_ec ce,
               ligand_to_prox_chain_id lpc
         WHERE r.ec = ce.ec
               AND lpc.ligand_id = lc.ligand_id AND lpc.chain_id = ce.chain_id
               AND lc.chebi_id = ANY(r.products)
       ) sq
 WHERE is_product = false AND sq.ligand_id = l.ligand_id;

 UPDATE credo_dev.ligands l
   SET is_cofactor = true
  FROM (
         SELECT DISTINCT lc.ligand_id
           FROM intenz.enzymes r, ligand_to_chebi_id lc, chain_id_to_ec ce,
                ligand_to_prox_chain_id lpc
          WHERE r.ec = ce.ec
                AND lpc.ligand_id = lc.ligand_id AND lpc.chain_id = ce.chain_id
                AND lc.chebi_id = ANY(r.cofactors)
       ) sq
 WHERE sq.ligand_id = l.ligand_id;
