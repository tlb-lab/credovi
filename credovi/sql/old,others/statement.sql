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
