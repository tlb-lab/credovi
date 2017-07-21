if [ "$2" ]
then
find /tlbnas/store/credo/data/ -mtime $2 -name ligands.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_ligands FROM \'{}\'\"
find /tlbnas/store/credo/data/ -mtime $2 -name aromaticrings.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_aromatic_rings FROM \'{}\'\"
find /tlbnas/store/credo/data/ -mtime $2 -name pisystems.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_pi_groups FROM \'{}\'\"
find /tlbnas/store/credo/data/ -mtime $2 -name atoms.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_atoms FROM \'{}\'\"
find /tlbnas/store/credo/data/ -mtime $2 -name contacts.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_contacts FROM \'{}\'\"
else
find /tlbnas/store/credo/data/ -name ligands.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_ligands FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name aromaticrings.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_aromatic_rings FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name pisystems.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_pi_groups FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name atoms.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_atoms FROM \'{}\'\"
find /tlbnas/store/credo/data/ -name contacts.credo | parallel psql cryst -p 5432 -c \"COPY $1.raw_contacts FROM \'{}\'\"
fi
