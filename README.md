## larda_cloud_sniffer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4723824.svg)](https://doi.org/10.5281/zenodo.4723824)

Identify coherent cloud objects in the Cloudnet classification and extract characteristic properties from layered mixed-phase clouds.


### sniffer_code

- the cloud sniffer is originally based on Bühl [2016 ACP]
- adapted for the new datasets, additinal measurements, new server and especially larda3

### cloud_properties

`.dat` files produced by `python3 cc_sniffer_ac.py --campaign lacros_dacapo --date 20181128` for a single day.
Contains all the `features` (range chunks of profiles) connected into clouds.

### cloud_collections 

Statistics of the single cloud features for the full campaign in a `.csv` table.
Created by `python3 cc_collector_ac.py --campaign lacros_dacapo --date 20181128`

The cloud collections for the LACROS campaigns at Leipzig, Limassol and Punta Arenas are included in the zenodo repository,
but **not in the github repository**.


### analysis_code

Collection of ipython notebooks to generate the analysis and the plots used in the recent publication.



### References
- Bühl, J., Seifert, P., Myagkov, A., and Ansmann, A.: Measuring ice- and liquid-water properties in mixed-phase cloud layers at the Leipzig Cloudnet station, Atmos. Chem. Phys., 16, 1060-10620, <https://doi.org/10.5194/acp-16-10609-2016>, 2016. 


### License
See the LICENSE file for more information
Copyright 2021, Martin Radenz
[MIT License](http://www.opensource.org/licenses/mit-license.php)


