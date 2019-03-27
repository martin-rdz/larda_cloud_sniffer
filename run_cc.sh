#!/bin/bash

camp='lacros_dacapo'

for i in 201812{01..31}; do 
#for i in 201901{01..31}; do 
#for i in 201902{01..28}; do 
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    #python3 cc_collector_ac.py --campaign lacros_dacapo --date $i --type mixed-phase
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
done

for i in 201901{01..31}; do 
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
done
for i in 201902{01..28}; do 
    python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
    python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
done


camp='lacros_cycare'
d="20161018"
while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "20180325" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign lacros_cycare --date $d;
    python3 cc_collector_ac.py --campaign lacros_cycare --date $d
done
