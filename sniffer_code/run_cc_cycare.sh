#!/bin/bash

#abort if a commadn fails
#set -e


#d="20170309"
#python3 cc_sniffer_ac.py --campaign lacros_cycare --date $d;
#python3 cc_collector_ac.py --campaign lacros_cycare --date $d

camp='lacros_cycare'
d="20161018"
#d="20161019"
#d="20161225"
#d="20170228"
#d="20170307"
#while [ "$(date -d "$d" +%Y%m%d)" -lt "$(date -d "20180325" +%Y%m%d)" ]; do
#    python3 cc_sniffer_ac.py --campaign lacros_cycare --date $d;
#    python3 cc_collector_ac.py --campaign lacros_cycare --date $d
#    d=$(date -d "$d + 1 day" +%Y%m%d)
#done

file="avail_dates_lacros_cycare.dat"

cat ${file} | while read i || [[ -n $i ]]; do
#while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20200630" +%Y%m%d)" ]; do
    echo ${i}
    if [[ "$i" != *"#"* ]]; then
        #python3 cc_sniffer_ac.py --campaign lacros_cycare --date $i;
        python3 cc_collector_ac.py --campaign lacros_cycare --date $i
    fi
    #i=$(date -d "$i + 1 day" +%Y%m%d)
done
exit 1
