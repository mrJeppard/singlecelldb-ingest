#!/bin/bash

FILES=/home/duncan/work/clusterdb-ingest/marker_jsons/*

SCR="python /home/duncan/work/clusterdb-ingest/ingest/marker_vals_gzipjson_in.py" 

for f in $FILES
do
  echo "Processing $f file..."
  $SCR -gz $f -db /home/duncan/work/clusterdb-ingest/cluster.db.ebi.filled
done
