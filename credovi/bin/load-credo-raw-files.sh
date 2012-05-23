find /merlin/Store/CREDO/data -name ligands.credo | parallel psql cryst -c \"COPY $1.raw_ligands FROM \'{}\'\"
find /merlin/Store/CREDO/data -name aromaticrings.credo | parallel psql cryst -c \"COPY $1.raw_aromatic_rings FROM \'{}\'\"
find /merlin/Store/CREDO/data -name atoms.credo | parallel psql cryst -c \"COPY $1.raw_atoms FROM \'{}\'\"
find /merlin/Store/CREDO/data -name contacts.credo | parallel psql cryst -c \"COPY $1.raw_contacts FROM \'{}\'\"
