-- CREATE A MAPPING BETWEEN LIGANDS AND FRAGMENTS
DO $$
    DECLARE
        lig_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL LIGANDS THAT DO NOT HAVE CLASHES AND ARE NOT DISORDERED
        FOR lig_id, biomol_id IN   SELECT ligand_id, biomolecule_id
                                     FROM credo.ligands
                                    WHERE is_disordered = false AND is_clashing = false
                                 ORDER BY 1
        LOOP
            EXECUTE
            '
            WITH ligand_fragments AS
                 (
                     INSERT INTO credo.ligand_fragments(biomolecule_id, ligand_id, ligand_component_id, fragment_id, hit)
                     SELECT DISTINCT r.biomolecule_id, lc.ligand_id, lc.ligand_component_id, fragment_id, hit
                       FROM credo.ligand_components lc
                       JOIN credo.residues r ON r.residue_id = lc.residue_id
                       JOIN pdbchem.chem_comp_fragments cf ON r.res_name = cf.het_id
                       JOIN pdbchem.chem_comp_fragment_atoms cfa ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                      WHERE lc.ligand_id = $1
                   ORDER BY 1,2,3
                  RETURNING credo.ligand_fragments.*
                 ),
                 ligand_fragment_atoms AS
                 (
                     INSERT INTO credo.ligand_fragment_atoms(ligand_fragment_id, atom_id)
                     SELECT DISTINCT lf.ligand_fragment_id, a.atom_id
                       FROM ligand_fragments lf
                            -- get the ligand component atoms
                       JOIN credo.ligand_components lc ON lc.ligand_component_id = lf.ligand_component_id
                       JOIN credo.residues r ON r.residue_id = lc.residue_id
                       JOIN credo.atoms a on a.residue_id = lc.residue_id
                            -- link them to the fragment atoms
                       JOIN pdbchem.chem_comp_fragments cf
                            ON cf.het_id = r.res_name AND cf.fragment_id = lf.fragment_id
                       JOIN pdbchem.chem_comp_fragment_atoms cfa
                            ON cfa.chem_comp_fragment_id = cf.chem_comp_fragment_id
                            AND lf.hit = cfa.hit
                            AND cfa.pdb_name = a.atom_name
                      WHERE lc.ligand_id = lf.ligand_id AND a.biomolecule_id = $2
                   ORDER BY 1,2
                  RETURNING credo.ligand_fragment_atoms.*
                 )
                 SELECT NULL;
            ' USING lig_id, biomol_id;
            RAISE NOTICE 'inserted ligand fragments and their atoms for ligand %', lig_id;
        END LOOP;
END$$;

CREATE  INDEX idx_ligand_fragments_biomolecule_id
        ON credo.ligand_fragments
        USING btree (biomolecule_id) WITH(FILLFACTOR=100);

CREATE  INDEX idx_ligand_fragments_ligand_id
        ON credo.ligand_fragments
        USING btree (ligand_id) WITH(FILLFACTOR=100);

CREATE  INDEX idx_ligand_fragments_ligand_component_id
        ON credo.ligand_fragments
        USING btree (ligand_component_id) WITH(FILLFACTOR=100);

CREATE  INDEX idx_ligand_fragments_fragment_id
        ON credo.ligand_fragments
        USING btree (fragment_id) WITH(FILLFACTOR=100);

CREATE  UNIQUE INDEX idx_ligand_fragment_atoms_ligand_fragment_id
        ON credo.ligand_fragment_atoms
        USING btree (ligand_fragment_id, atom_id) WITH(FILLFACTOR=100);

CREATE  INDEX idx_ligand_fragment_atoms_atom_id
        ON credo.ligand_fragment_atoms
        USING btree (atom_id) WITH(FILLFACTOR=100);


-- SET A FLAG OF ROOT LIGAND FRAGMENTS
UPDATE credo.ligand_fragments lf
   SET is_root = true
  FROM credo.ligand_components lc,
       pdbchem.fragment_hierarchies fh
 WHERE lc.ligand_component_id = lf.ligand_component_id
       AND fh.parent_id = lf.fragment_id
       AND fh.order_parent = 0
       AND fh.het_id = lc.het_id;

DO $$
    DECLARE
        lig_frag_id INTEGER;
        biomol_id INTEGER;
    BEGIN
        -- LOOP THROUGH ALL BIOMOLECULES
        FOR lig_frag_id, biomol_id IN   SELECT lf.ligand_fragment_id, lf.biomolecule_id
                                          FROM credo.ligand_fragments lf
                                          JOIN credo.ligands USING(ligand_id)
                                      ORDER BY 1
        LOOP
            EXECUTE
                '
                   WITH fragment_contacts AS
                        (
                           SELECT fc.ligand_fragment_id,
                                  COUNT(DISTINCT fc.atom_id) AS num_int_atoms,
                                  SUM(fc.is_covalent::int) as num_covalent,
                                  SUM(fc.is_vdw_clash::int) as num_vdw_clash,
                                  SUM(fc.is_vdw::int) AS num_vdw,
                                  SUM(fc.is_proximal::int) AS num_proximal,
                                  SUM(fc.is_hbond::int) AS num_hbond,
                                  SUM(fc.is_weak_hbond::int) AS num_weak_hbond,
                                  SUM(fc.is_xbond::int) as num_xbond,
                                  SUM(fc.is_ionic::int) as num_ionic,
                                  SUM(fc.is_metal_complex::int) AS num_metal_complex,
                                  SUM(fc.is_aromatic::int) AS num_aromatic,
                                  SUM(fc.is_hydrophobic::int) AS num_hydrophobic,
                                  SUM(fc.is_carbonyl::int) AS num_carbonyl
                             FROM (
                                      SELECT lfa.ligand_fragment_id,
                                             cs.atom_bgn_id as atom_fragment_id, cs.atom_end_id as atom_id, cs.*
                                        FROM credo.ligand_fragment_atoms lfa
                                             -- FRAGMENT ATOM IS BGN
                                        JOIN credo.contacts cs ON cs.atom_bgn_id = lfa.atom_id
                                       WHERE lfa.ligand_fragment_id = $1
                                             AND cs.biomolecule_id = $2
                                             AND cs.is_same_entity = false
                                             AND cs.distance <= 4.5
                                             -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                             AND cs.structural_interaction_type_bm & 56 > 0
                                   UNION ALL
                                      SELECT lfa.ligand_fragment_id,
                                             cs.atom_end_id as atom_fragment_id, cs.atom_bgn_id as atom_id, cs.*
                                        FROM credo.ligand_fragment_atoms lfa
                                             -- FRAGMENT ATOM IS END
                                        JOIN credo.contacts cs ON cs.atom_end_id = lfa.atom_id
                                       WHERE lfa.ligand_fragment_id = $1
                                             AND cs.biomolecule_id = $2
                                             AND cs.is_same_entity = false
                                             AND cs.distance <= 4.5
                                             -- INTERACTION MUST BE WITH POLYMER ATOM (PROT/DNA/RNA)
                                             AND cs.structural_interaction_type_bm & 3584 > 0
                                  ) fc
                         GROUP BY fc.ligand_fragment_id
                         ORDER BY fc.ligand_fragment_id
                        )
                 UPDATE credo.ligand_fragments lf
                    SET num_int_atoms = fc.num_int_atoms,
                        num_contacts = fc.num_covalent + fc.num_vdw_clash + fc.num_vdw + fc.num_proximal,
                        num_covalent = fc.num_covalent,
                        num_vdw_clash = fc.num_vdw_clash,
                        num_vdw = fc.num_vdw,
                        num_proximal = fc.num_proximal,
                        num_hbond = fc.num_hbond,
                        num_weak_hbond = fc.num_weak_hbond,
                        num_xbond = fc.num_xbond,
                        num_ionic = fc.num_ionic,
                        num_metal_complex = fc.num_metal_complex,
                        num_aromatic = fc.num_aromatic,
                        num_hydrophobic = fc.num_hydrophobic,
                        num_carbonyl = fc.num_carbonyl
                   FROM fragment_contacts fc
                  WHERE fc.ligand_fragment_id = lf.ligand_fragment_id;
                  ' USING lig_frag_id, biomol_id;

                RAISE NOTICE 'updated fragment contact densities for ligand fragment %', lig_frag_id;
        END LOOP;
END$$;

UPDATE credo.ligand_fragments u
   SET fcd = sq.fcd, fad = sq.fad
  FROM (
        -- get the number of ligand heavy atoms and interacting atoms from the root fragments
        WITH ligands AS
        (
         SELECT ligand_id,
                num_hvy_atoms, num_int_atoms, num_contacts,
                num_int_atoms / num_hvy_atoms::numeric AS int_ratio,
                num_contacts / num_hvy_atoms::numeric AS cs_ratio
           FROM credo.ligand_fragments lf
           JOIN credo.ligands l USING(ligand_id)
          WHERE lf.is_root = true
        )
        SELECT lf.ligand_fragment_id,
               (lf.num_int_atoms / f.num_hvy_atoms::numeric) / l.int_ratio as fad,
               (lf.num_contacts / f.num_hvy_atoms::numeric) / l.cs_ratio as fcd
          FROM credo.ligand_fragments lf
          JOIN ligands l ON l.ligand_id = lf.ligand_id
          JOIN pdbchem.fragments f ON f.fragment_id = lf.fragment_id
       ) sq
 WHERE u.ligand_fragment_id = sq.ligand_fragment_id;
