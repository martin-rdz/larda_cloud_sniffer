#!/bin/bash

#abort if a commadn fails
#set -e



camp='davos'
d="20190210"
e="20190321"
while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "$e" +%Y%m%d)" ]; do
    python3 cc_sniffer_ac.py --campaign davos --date $d;
    python3 cc_collector_ac.py --campaign davos --date $d
    d=$(date -d "$d + 1 day" +%Y%m%d)
done
