DROP TABLE IF EXISTS credo.lig_nuc_interactions;
CREATE TABLE credo.lig_nuc_interactions (
    lig_nuc_interaction_id serial NOT NULL,
    biomolecule_id integer NOT NULL,
    ligand_id integer NOT NULL,
    chain_nuc_id integer NOT NULL,
    path ptree,
    is_quaternary boolean NOT NULL DEFAULT false,
    CONSTRAINT lig_nuc_interactions_pkey PRIMARY KEY (lig_nuc_interaction_id)
);

CREATE INDEX idx_lig_nuc_interactions_biomolecule_id
       ON credo.lig_nuc_interactions USING btree (biomolecule_id);

CREATE INDEX idx_lig_nuc_interactions_ligand_id
       ON credo.lig_nuc_interactions USING btree (ligand_id);

CREATE INDEX idx_lig_nuc_interactions_chain_nuc_id
       ON credo.lig_nuc_interactions USING btree (chain_nuc_id);

-- insert new ligand-nucleic acid interactions
DO $$
    DECLARE
        biomol_id INTEGER;
    BEGIN
        -- loop through all biomolecules that are not in the ligand-nucleic acid interactions table yet
        FOR biomol_id IN   SELECT DISTINCT biomolecule_id
                             FROM credo.chains
                             JOIN credo.oligonucleotides USING(chain_id)
                            WHERE biomolecule_id > (SELECT COALESCE(max(biomolecule_id),0)
                                                      FROM credo.lig_nuc_interactions)
                         ORDER BY 1
        LOOP
            EXECUTE
            '
              INSERT INTO credo.lig_nuc_interactions(biomolecule_id, ligand_id, chain_nuc_id)
              SELECT DISTINCT rip.biomolecule_id, lc.ligand_id, n.chain_id
                FROM credo.residue_interaction_pairs rip
                JOIN credo.nucleotides n ON n.residue_id = rip.residue_bgn_id
                JOIN credo.ligand_components lc ON lc.residue_id = rip.residue_end_id
               WHERE rip.biomolecule_id = $1
               UNION
              SELECT DISTINCT rip.biomolecule_id, lc.ligand_id, n.chain_id
                FROM credo.residue_interaction_pairs rip
                JOIN credo.nucleotides n ON n.residue_id = rip.residue_end_id
                JOIN credo.ligand_components lc ON lc.residue_id = rip.residue_bgn_id
               WHERE rip.biomolecule_id = $1
            ORDER BY 1,2,3;
            ' USING biomol_id;

            RAISE NOTICE 'inserted ligand-nucleic interactions for biomolecule %', biomol_id;
        END LOOP;
END$$;

--
UPDATE credo.lig_nuc_interactions li
   SET path = sq.path
  FROM (
        SELECT li.lig_nuc_interaction_id,
                (
                 b.path || ('LNI:' || li.lig_nuc_interaction_id) || l.ligand_name || c.pdb_chain_id
                ) ::ptree AS path
          FROM credo.lig_nuc_interactions li
          JOIN credo.biomolecules b USING(biomolecule_id)
          JOIN credo.ligands l ON li.ligand_id = l.ligand_id
          JOIN credo.chains c ON li.chain_nuc_id = c.chain_id
       ) sq
 WHERE sq.lig_nuc_interaction_id = li.lig_nuc_interaction_id;