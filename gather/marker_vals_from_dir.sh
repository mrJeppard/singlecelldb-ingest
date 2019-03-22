#!/bin/bash

FILES=/home/duncan/work/ebiscanpy-clusterdb/data/*
OUTDIR=/home/duncan/work/clusterdb-ingest/marker_jsons

SCR="python3 /home/duncan/work/clusterdb-ingest/gather/marker_genes/markervalues.py"
ENV=/home/duncan/work/clusterdb-ingest/env

source $ENV/bin/activate

for f in $FILES
do
  echo "Processing $f file..."
  $SCR -ia $f -cs sc3_cluster_solutions -o $OUTDIR
done