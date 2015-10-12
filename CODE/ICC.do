* upload either the household level data or individual level data
use "DATA\q1q4_ctrl.dta", clear
use "Data\q2_ctrl.dta", clear
* format
destring cluster, replace
destring dataid, replace
* get ICC for each exposure
gllamm q11_tube, family(binomial) link(logit) i(cluster) eform adapt 
nlcom [clus1]_cons^2 / ( [clus1]_cons^2 +3.14*3.14/3)



