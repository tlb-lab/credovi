for path in /tlbnas/mirror/ensembl/homo_sapiens_core_70_37/*.txt.gz
do
    filename=`basename $path`
    table=${filename%%.*}
    filesize=$(stat -c%s "$path")

    if [[ $filesize -lt 5242900000 ]]; then
    	echo "loading ${filename} into table ${table}..."

    	# COPY TABLE FROM STDIN
    	zcat $path | sed 's/0000-00-00 00:00:00/\\N/g' | psql -p 5432 -d cryst -c "COPY ensemblcore.${table} FROM STDIN"
    else
       echo "${filename} is too large (${filesize}) and will be ignored." 
    fi
done
