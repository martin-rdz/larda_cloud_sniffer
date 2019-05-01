#!/bin/bash

#abort if a command fails
#set -e


#camp='lacros_dacapo'
#i='20181204'
i='20181211'
python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
exit 1

#for i in 201812{01..31}; do 
#for i in 201901{01..31}; do 
#for i in 201902{01..28}; do 
#    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    #python3 cc_collector_ac.py --campaign lacros_dacapo --date $i --type mixed-phase
#    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
#done

i="20181128"
#skip the 25jan
#i="20190201"
while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20190415" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
    i=$(date -d "$i + 1 day" +%Y%m%d)
done
exit 1

camp='lacros_cycare'
d="20161018"
#d="20161019"
while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "20180325" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign lacros_cycare --date $d;
    python3 cc_collector_ac.py --campaign lacros_cycare --date $d
    d=$(date -d "$d + 1 day" +%Y%m%d)
done
