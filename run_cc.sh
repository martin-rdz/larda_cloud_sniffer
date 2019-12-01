#!/bin/bash

#abort if a command fails
#set -e


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
#i="20191028"
#skip the 25jan
#i="20190201"
while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20191104" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
    i=$(date -d "$i + 1 day" +%Y%m%d)
done
exit 1

