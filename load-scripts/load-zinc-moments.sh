#!/bin/bash
FILES=/merlin/Store/ZINC/usrcat/*.usrcat.gz

for f in $FILES
do
  echo "Processing $f file..."
  zcat $f | psql -p 5433 -c "COPY zinc._moments FROM STDIN" -d cryst
  psql -p 5433 -d cryst -c "INSERT INTO zinc.moments SELECT * FROM zinc._moments"
  psql -p 5433 -d cryst -c "TRUNCATE zinc._moments"
  
done
