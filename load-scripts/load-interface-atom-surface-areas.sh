find /merlin/Store/CREDO/data -name interface_atom_surface_areas.credo | parallel psql cryst -c \"COPY credo.raw_interface_atom_surface_areas FROM \'{}\'\"
