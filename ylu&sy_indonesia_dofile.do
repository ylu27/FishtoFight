**********Replication data for WD-25724**********
*packages
ssc install ivreg2
ssc install spatial_hac_iv
net install st0292
//Please also find attached .ado files in respiratory

*Table 2. Summary statistics
use "stata\ylu&sy_indonesia_sample.dta"  

summarize conflict conflict50nmi conflict100nmi conflict150nmi conflict200nmi fatalities catch000 indcatch000 nidcatch000 iuu000 ind_iuu000 nid_iuu000 ocean_prd

clear all
//------------------------



*Table 3-A. First-stage coefficients and Conley-HAC standard error 

use "stata\ylu&sy_indonesia_sample.dta"  
//distance cut-off at 111 km (determined by 1x1-degree cell size); lagcutoff at 1 year because of cross-sectional analysis (our study has no time lag)
//same setting for the rest estimates
ols_spatial_HAC indcatch000 ocean_prd prodummy*, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar
ols_spatial_HAC nidcatch000 ocean_prd prodummy*, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar
ols_spatial_HAC catch000 ocean_prd prodummy*, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar

*Table 3-B. Diagnostics (please find Kleibergen-Paap LM statistics and F-statistics of excluded instruments from first-stage panel)
ivreg2 conflict100nmi (indcatch000 = ocean_prd) prodummy*, robust first noconstant
ivreg2 conflict100nmi (nidcatch000 = ocean_prd) prodummy*, robust first noconstant
ivreg2 conflict100nmi (catch000 = ocean_prd) prodummy*, robust first noconstant

clear all
//-----------------------



*Table 4. Effects of industrial and non-industrial fisheries on the number of conflicts and fatalities within 100 nautical miles
use "stata\ylu&sy_indonesia_sample.dta"  

//Perliminary results - pooled OLS without province fixed effects
ols_spatial_HAC conflict100nmi indcatch000, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar
ols_spatial_HAC conflict100nmi nidcatch000, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar
ols_spatial_HAC conflict100nmi catch000, p(OBJECTID) t(year) lat(latitude) lon(longitude) distcutoff(111) lagcutoff(1) dropvar

//Core specification (Table 4) - pooled 2sls without province fixed effects
drop if missing(ocean_prd) //drop missing value in ocean_productivity

spatial_hac_iv conflict100nmi (indcatch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi (nidcatch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi (catch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

//Core results (Table 4 cont.) - 2sls with province fixed effects
spatial_hac_iv conflict100nmi prodummy* (indcatch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* (nidcatch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* (catch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv fatalities100nmi prodummy* (catch000 = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 5. Spatial autoregressive regressions
use "stata\ylu&sy_indonesia_sample.dta"  

//Create spectral-normalized spatial weight matrix based on the haversine distance (in miles) for the longitude and latitude of sample cells
drop if missing(ocean_prd)
spmat idistance obj1 longitude latitude, id(OBJECTID) dfunction(dhaversine, miles) normalize(spec) 

//Re-estimation of baseline results in Table 4 
//elmat ~ extent of spatial interactions in the dependent variable 
//dlmat ~ extent of spatial interactions in the error-term
spivreg conflict100nmi prodummy* (indcatch000 = ocean_prd), id(OBJECTID)  elmat(obj1) dlmat(obj1) noconstant hetero
spivreg conflict100nmi prodummy* (nidcatch000 = ocean_prd), id(OBJECTID)  elmat(obj1) dlmat(obj1) noconstant hetero
spivreg conflict100nmi prodummy* (catch000 = ocean_prd), id(OBJECTID)  elmat(obj1) dlmat(obj1) noconstant hetero
spivreg fatalities100nmi prodummy* (catch000 = ocean_prd), id(OBJECTID) elmat(obj1) dlmat(obj1) noconstant hetero

clear all
//-----------------------



*Table 6. Lagged effects of fisheries on the number of conflicts 
use "stata\ylu&sy_indonesia_sample.dta"  

//Variable definition: annual_catch1014 ~ mean quantity of fish caught annually during 2010 to 2014; ocp1014mean ~ mean ocean productivity between 2010 and 2014
drop if missing(annual_catch1014)
drop if missing(ocp1014mean)
spatial_hac_iv conflict100nmi prodummy* (annual_indcatch1014 = ocp1014mean), robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* (annual_nidcatch1014 = ocp1014mean), robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* (annual_catch1014 = ocp1014mean),  robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 7. Lagged effects of fisheries on the number of conflicts within the year
use "stata\ylu&sy_indonesia_sample.dta"  

//Variable definition: conflictm712 ~ number of conflicts from July to December in 2015; ocp16 ~ ocean productivity from January to June in 2015 
drop if missing(ocp16)
spatial_hac_iv conflictm712 prodummy* (indcatch000 = ocp16), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflictm712 prodummy* (nidcatch000 = ocp16),  robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflictm712 prodummy* (catch000 = ocp16), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 8. Persistence in conflict
use "stata\ylu&sy_indonesia_sample.dta"  

//Variable definition: conflict2016 ~ number of conflicts occurred in 2016
drop if missing(ocean_prd)
spatial_hac_iv conflict2016 (indcatch000 = ocean_prd) conflict100nmi prodummy*, robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict2016 ( nidcatch000 = ocean_prd) conflict100nmi prodummy*, robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict2016 (catch000 = ocean_prd) conflict100nmi prodummy*, robust lon(longitude) lat(latitude) timevar(year) panelvar(OBJECTID) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 9. Effects of catch variability anomalies on conflict
use "stata\ylu&sy_indonesia_sample.dta"  

//Variable definition: catchdev ~ catch variability anomalies
drop if missing(ocean_prd)
spatial_hac_iv conflict100nmi prodummy* ( indcatchdev = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( nidcatchdev = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( catchdev = ocean_prd), lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 10. Catch landings and conflict by conflict types
use "stata\ylu&sy_indonesia_sample.dta"  

//Effect of fishery catches on conflicts at different level of violence (i.e., protest; civilian; battle)
drop if missing(ocean_prd)
spatial_hac_iv protest_count100nmi prodummy* (catch000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar 
spatial_hac_iv civilian_count100nmi prodummy* (catch000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar 
spatial_hac_iv battle_count100nmi prodummy* (catch000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar 

clear all
//-----------------------



*Table 11. Regional differences in the impact of fisheries on the number of conflicts 
use "stata\ylu&sy_indonesia_sample.dta"  

//We use devgroup 1,2,3,4 (i.e., A,B,C,D) to run regressions with sub-region-sample
drop if missing(ocean_prd)
spatial_hac_iv conflict100nmi prodummy* ( catch000 = ocean_prd) if devgroup1==1 , robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( catch000 = ocean_prd) if devgroup2==1 , robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( catch000 = ocean_prd) if devgroup3==1 , robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( catch000 = ocean_prd) if devgroup4==1 , robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Table 12. Impact of IUU fishing
use "stata\ylu&sy_indonesia_sample.dta"  

drop if missing(ocean_prd)
spatial_hac_iv conflict100nmi prodummy* (ind_iuu000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( nid_iuu000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar
spatial_hac_iv conflict100nmi prodummy* ( iuu000 = ocean_prd), robust lon(longitude) lat(latitude) timevar(year) panelvar( OBJECTID ) lagcutoff(1) distcutoff(111) dropvar

clear all
//-----------------------



*Appendix

*Table C1. Correlation of variables
use "stata\ylu&sy_indonesia_sample.dta"  

corr conflict conflict50nmi conflict100nmi conflict150nmi conflict200nmi fatalities catch000 indcatch000 nidcatch000 iuu000 ind_iuu000 nid_iuu000 ocean_prd

clear all
//-----------------------



*Table E. Poisson IV with a control function estimator
use "stata\ylu&sy_indonesia_sample.dta"  

//Re-estimation of baseline results by Possion-IV with heterosckedasiticity-robust standard error
ivpoisson cfunction conflict100nmi prodummy* (indcatch000 = ocean_prd), vce(robust)
ivpoisson cfunction conflict100nmi prodummy* (nidcatch000 = ocean_prd), vce(robust)
ivpoisson cfunction conflict100nmi prodummy* (catch000 = ocean_prd), vce(robust)

clear all


**********End of Replication**********
