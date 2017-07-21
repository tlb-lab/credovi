CREATE OR REPLACE FUNCTION "credo_dev"."entity_bm_to_string"("entity_type_bm" int4) 
    RETURNS "pg_catalog"."text" AS 
$$
DECLARE
    entt_bm_tmp int4;
-- Entity_type_bm: solvent=1, ligand=2, saccharide=4, rna=8, dna=16, protein=32
    entt_names text[] := array['SOLV','LIG','SUG','RNA','DNA','PROT'];
    entt_type_arr text[];
BEGIN
  FOR i IN 1 .. 6 LOOP
    entt_bm_tmp := 1 << i-1;
    --RAISE DEBUG '% => % (%)', i, entt_bm_tmp, entt_names[i];
    IF entity_type_bm & entt_bm_tmp = entt_bm_tmp THEN
        entt_type_arr := entt_type_arr || entt_names[i];
    END IF;
  END LOOP;
  RETURN array_to_string(entt_type_arr, '|');
END
$$
LANGUAGE 'plpgsql' IMMUTABLE;


CREATE OR REPLACE FUNCTION "credo_dev"."struct_inter_bm_to_string"(structural_interaction_type_bm "int4")
  RETURNS "pg_catalog"."text" AS $$
DECLARE
	ntt_bm_bgn int4;
	ntt_bm_end int4;
	ntt_bgn_str text;
	ntt_end_str text;
BEGIN
  -- Entity_type_bm: solvent=1, ligand=2, saccharide=4, rna=8, dna=16, protein=32        
  ntt_bm_bgn := structural_interaction_type_bm >> 6;
	ntt_bm_end := structural_interaction_type_bm & 63; -- B'000000111111'::int4

	IF ntt_bm_bgn = 0 OR ntt_bm_end = 0 THEN
		RAISE EXCEPTION '% is not a valid structural interaction type bitmask', structural_interaction_type_bm;
	ELSE
		ntt_bgn_str := entity_bm_to_string(ntt_bm_bgn);
		ntt_end_str := entity_bm_to_string(ntt_bm_end);
		-- RAISE NOTICE '%s - BM1 (%): %, BM2 (%): %', structural_interaction_type_bm::bit(12), ntt_bm_bgn, ntt_bgn_str, ntt_bm_end, ntt_end_str;
	END IF;
	RETURN ntt_bgn_str || '-' || ntt_end_str;
END
$$
LANGUAGE 'plpgsql' IMMUTABLE;