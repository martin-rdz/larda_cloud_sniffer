## larda_cloud_sniffer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4723824.svg)](https://doi.org/10.5281/zenodo.4723824)

Identify coherent cloud objects in the Cloudnet classification and extract characteristic properties from layered mixed-phase clouds.

The structure of the folders is as following:

### sniffer_code

- the cloud sniffer is originally based on Bühl [2016 ACP]
- adapted for the new datasets, additinal measurements, new server and especially [larda3](https://github.com/lacros-tropos/larda).

### cloud_properties

`.dat` files produced by `python3 cc_sniffer_ac.py --campaign lacros_dacapo --date 20181128` for a single day.
Contains all the `features` (range chunks of profiles) connected into clouds.

### cloud_collections 

Statistics of the single cloud features for the full campaign in a `.csv` table.
Created by `python3 cc_collector_ac.py --campaign lacros_dacapo --date 20181128`

The cloud collections for the LACROS campaigns at Leipzig, Limassol and Punta Arenas are included in the zenodo repository,
but **not in the github repository**.
To obtain those `git clone` this repository and then get the `.csv` files from the [zenodo repository](https://doi.org/10.5281/zenodo.4723823).


### analysis_code

Collection of ipython notebooks to generate the analysis and the plots used in the recent publication.



### References
- Bühl, J., Seifert, P., Myagkov, A., and Ansmann, A.: Measuring ice- and liquid-water properties in mixed-phase cloud layers at the Leipzig Cloudnet station, Atmos. Chem. Phys., 16, 1060-10620, <https://doi.org/10.5194/acp-16-10609-2016>, 2016. 
- Radenz, M., Bühl, J., Seifert, P., Baars, H., Engelmann, R., Barja González, B., Mamouri, R.-E., Zamorano, F., and Ansmann, A.: Hemispheric contrasts in ice formation in stratiform mixed-phase clouds: Disentangling the role of aerosol and dynamics with ground-based remote sensing, Atmos. Chem. Phys. Discuss. [accepted for publication], <https://doi.org/10.5194/acp-2021-360>, 2021. 


### License
See the LICENSE file for more information
Copyright 2021, Martin Radenz
[MIT License](http://www.opensource.org/licenses/mit-license.php)


