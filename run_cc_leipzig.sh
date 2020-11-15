#!/bin/bash

#abort if a command fails
set -e


i="20150708"
# skip 20150805 - 0810
i="20150811"
# skip 20150919 - 0920
i="20150921"
# skip 20151001 - 1013
i="20151014"
# skip 20151031 - 1101
i="20151102"
e="20180904"
# skip 20150330 - 0427
i="20160428"
e="20160701"

i="20160702"
e="20160720"

i="20160726"
e="20160730"

i="20160808"
e="20160813"

i="20160816"
e="20160820"

i="20160823"
e="20160923"

i="20180428"
e="20180718"

i="20180719"
e="20180730"

i="20180803"
e="20180807"


i="20150708"
e="20180807"


i="20140101"
e="20180807"

cat avail_dates_lacros_leipzig.dat | while read i || [[ -n $i ]]; do
#while [ "$(date -d "$i" +%Y%m%d)" -lt "$(date -d "$e" +%Y%m%d)" ]; do
    if [[ "$i" != *"#"* ]]; then
       	python3 cc_sniffer_ac.py --campaign lacros_leipzig --date $i;
        python3 cc_collector_ac.py --campaign lacros_leipzig --date $i
    fi
    #i=$(date -d "$i + 1 day" +%Y%m%d)
done
