find /tlbnas/store/credo/data/ -name binding_site_atom_surface_areas.credo | parallel psql cryst -c \"COPY $1.raw_binding_site_atom_surface_areas FROM \'{}\'\"
