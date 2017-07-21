=========================================================
CREDOVI: Software to create and manage the CREDO database
=========================================================


CREDO Creation/Update Instructions

==== Broad Overview ====

Phase 0 - Updating PDB source data (location: spunky, tlbnas):
    TLBNAS keeps an up-to-date MMCIF mirror on /tlbnas/mirror/pdb/data/structures/all/mmCIF. This needs to be loaded onto MySQL using the 'mmcif.rake' script on Spunky (TODO: cronjob).
    Sync the mmCIF mirror which is kept up to date on ''spunky'' via a cron job, with the Postgres database on ''bahamut'', using ''my2pgsql.py'' (skip atom_site and atom_site_anisotrop tables - not needed).	(Warning: The 'software' table contains 'NULL' values on the 'name' column which get converted to NULLs, violating the table constraint. Edit to NONE (or whatever) first.)
    TLBNAS should update /tlbnas/store/pdbe/pdbechem/ and /tlbnas/mirror/pdbe/sifts/ regularly, but make sure they are recent enough.
    
Phase 1 - Generating base and support schemas (location: any with OpenEye and access to DB and TLBNAS mirror, e.g. marid)
    - Extract the biomolecule transformations from PDB headers and store in the database, using ''create_biomt_from_header.py''. These data go into ''pdb_dev.biomt'' and ''pdb_dev.biomt_ops''
    - Create rest of 'pdb_dev' schema tables (pdb_dev.res_map and pdb_dev.map_regions, etc) with "create-res-map.py". Table 'banned' is not automatically generated and should be copied over from the previous version.
    NOTE: Since the last Adrian update, the EBI have changed their sifts XML format, which moved the NCBI annotation from MapRegion to Residue and apparently forgot to include the EC one anywhere altogether. This has consequences on CREDO's cross-reference tables, which rely on the information being on the map_region tables. 
	UPDATE 2016: EBI seems to have addressed this. However, entries between this and the official CSV SIFTS don't quite agree...
	- Create 'drugbank_dev' schema and populate using 'drugbank_loader.py' with the latest XML from Drugbank: http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip
		UPDATE 2017: The Drugbank XML can't be directly downloaded anymore, as it requires a login at the website, so the download needs to be done manually
	- Update 'emolecules' schema with 'emolecules_loader.py with the latest SMILES data from: http://downloads.emolecules.com/free/
    - Create 'pdbchem_dev' schema: 1) create_tables_from_cif.py > (copy astex_solvents table) > create_chem_comps.py > create-fragments.py  
	                               2) generate-conformers.py (see source code for parallel command) > create-chem-conformers.py (make sure usrcat package is installed)
	- Update 'sifts' schema from data on ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/. 
	  The Rake file 'msd_sifts.rake' takes care of this (requires Sequel package), but the CSVs can also be imported manually (although column names need tweaking)
    
Phase 2 - Generating CREDO's raw tables (location: multi-core server with OpenEye and access to DB and TLBNAS, e.g. marid):
        "python credovi.py mmcif currentpdbs | parallel --eta --halt 1 -n 6 python credovi.py credo preparepdb -Q --clean --oeb --pdb -S {1},{2},{3},{4},{5},{6}"
			(use -I instead of --clean to do an incremental update instead of generating the OEBs from scratch)
        "python credovi.py mmcif currentpdbs | parallel --eta --halt 1 -n 6 python credovi.py credo contacts -Q -S {1},{2},{3},{4},{5},{6}"  (should take about 8-12 hours)
			(use -I to do an incremental update instead of generating the contacts from scratch)
        "python credovi.py mmcif currentpdbs | parallel --eta --halt 1 -n 6 python credovi.py ligand surfareas -Q -S {1},{2},{3},{4},{5},{6}"
    
	Create tables on the schema with ''python credovi.py db create --sure --echo --checkfirst''
    
	
Phase 3 - Loading and generating CREDO's derived tables (SQL based)
	Run "bash load-credo-db-data.sh credo_dev" and "bash load-bindingsite-atom-surface-areas.sh credo_dev" on the database server (Bahamut)
    Execute ''populate.sql'' (e.g. 'psql -d cryst -f populate.sql -e' from Bahamut) (this will take a while)
	
	After the initial populating of the ligand related schemas, the following can be performed:    
		"python credovi.py mmcif currentpdbs | parallel --eta --halt 1 -n 6 python credovi.py ligand molstrings --usr -Q -S {1},{2},{3},{4},{5},{6}"
		Execute: "CLUSTER credo_dev.ligand_usr" on SQL after finishing the previous operation.

    Run the rest of the SQL scripts required, e.g. in PgAdmin:
            1) update-mmcif-dependent-fields.sql -> xrefs.sql -> update-ligands-with-kegg.sql, update-ligands-with-intenz.sql
             (intenz requires up-to-date pdb_chain_* from SIFTs schema, themselves created from data at ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/)
            2) ligand-fragments.sql, ligand-ligand-interactions.sql, ligand-nucleic-acid-interactions.sql
                        
    Copy from ''credo_dev'' to ''credo'', to overwrite the old CREDO

==== General Notes ====

    Keep an eye out for any schema changes in dependencies, e.g. ChEMBL
    The API and web interface should continue to work as long as the database schemas are not changed
    Keep an eye on the number of ''biomolecule_id''s. We'll need to increase ''current_biomol_max'' as the number increases to keep up with the number of partitions required for contacts tables.

