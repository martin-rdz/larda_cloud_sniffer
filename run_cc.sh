#!/bin/bash

#abort if a command fails
set -e
set -o pipefail


#camp='lacros_dacapo'
#i='20181204'
#i='20181211'
#python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
#python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
#exit 1

#for i in 201812{01..31}; do 
#for i in 201901{01..31}; do 
#for i in 201902{01..28}; do 
#    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    #python3 cc_collector_ac.py --campaign lacros_dacapo --date $i --type mixed-phase
#    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
#done

i="20181128"
i="20190111"
i="20190126"
i="20190327"
i="20190501"
i="20190605"
i="20190715"
i="20190923"
i="20191207"
i="20191215"
i="20200105"
i="20200123"
i="20200125"
i="20200128"
i="20200205"
i="20200212"
i="20200218"
i="20200304"
i="20200309"
i="20200312"
#i="20191028"
#skip the 25jan
#i="20190201"
cat avail_dates_lacros_dacapo.dat | while read i || [[ -n $i ]]; do
#while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20200630" +%Y%m%d)" ]; do
    echo ${i}
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
    #i=$(date -d "$i + 1 day" +%Y%m%d)
done
exit 1

