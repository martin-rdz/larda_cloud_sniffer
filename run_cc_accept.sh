#!/bin/bash

#abort if a commadn fails
#set -e


#d="20170309"
#python3 cc_sniffer_ac.py --campaign lacros_cycare --date $d;
#python3 cc_collector_ac.py --campaign lacros_cycare --date $d

camp='lacros_accept'
d="20141102"
while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "20141103" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign lacros_accept --date $d;
    python3 cc_collector_ac.py --campaign lacros_accept --date $d
    d=$(date -d "$d + 1 day" +%Y%m%d)
done
