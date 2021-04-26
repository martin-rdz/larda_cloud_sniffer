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
# 20190520 causes issues
i="20190522"
i="20190601"
# 20190714 causes issues
i="20190715"
# 20191105 causes issues
# the missing november 2019
i="20191206"
# 20191226 causes issues
i="20200105"
# 20200122 causes issues
i="20200123"
# 20200302 causes issues
i="20200304"
# 20200323 causes issues
i="20200331"
# 20200426 causes issues
# the missing may 2020
i="20200603"
# 20200622 causes issues
i="20200629"

# 20200704 causes issues
i="20200706"
# 20200722 causes issues
i="20200724"
# 20200822 causes issues
i="20200906"
# 20200924 causes issues
i="20200930"
# 20201124 causes issues
i="20201214"
# 20201228 strictly increasing error
i="20201229"


i="20200525"


file="avail_dates_lacros_dacapo_interim.dat"
file="avail_dates_lacros_dacapo_new.dat"
file="avail_dates_lacros_dacapo_new_inclMay.dat"

cat ${file}  | while read i || [[ -n $i ]]; do
#while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20210101" +%Y%m%d)" ]; do
#while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "20200603" +%Y%m%d)" ]; do
    echo ${i}
    if [[ "$i" != *"#"* ]]; then
        #python3 cc_sniffer_ac.py --campaign lacros_dacapo --date $i;
        python3 cc_collector_ac.py --campaign lacros_dacapo --date $i
    fi
    #i=$(date -d "$i + 1 day" +%Y%m%d)
done
exit 1

