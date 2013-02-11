find /tlbnas/store/credo/data/ -name ligands.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_ligands FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name aromaticrings.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_aromatic_rings FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name atoms.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_atoms FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name contacts.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_contacts FROM \'{}\'\"