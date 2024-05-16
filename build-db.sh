#!/bin/sh

DB_PATH="db/test"
build_db -k 35 -l 31 -S 1111111111111111101010101010101 \
         -H $DB_PATH/hash.k2d.tmp -t $DB_PATH/taxo.k2d.tmp -o $DB_PATH/opts.k2d.tmp \
		 -n $DB_PATH/taxonomy/ -m $DB_PATH/seqid2taxid.map -c 565654308 -p 1 \
		 -B 16384 -b 16384 -r 0