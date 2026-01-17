clear all
set matsize 11000
set maxvar 32000
set more off
set type double
cd "C:\Users\arpit\OneDrive - University of Pittsburgh\Pitt PhD Study Material\11 RESEARCH\Second Year\Marla Solving in Levels\Code\Code for SYP"



*************** INPUT CONCORDANCES AND CODES ***************************

import delimited using Numerics/Input_Files/iso_codes.csv, clear
rename v1 country_i
rename v2 iso_i
save Numerics/iso_codes_i, replace
rename country_i country_j
rename iso_i iso_j
save Numerics/iso_codes_j, replace


*** WIOD categories ***
import delimited using Numerics/Input_Files/Cat_WIOD_new.csv, clear
rename rnr code
rename industrycode rnr
rename industrydescription industry
keep industry code sector rnr
save Numerics/SectorsWIOD, replace


***********Get GO as weights to construct RoW values

use Numerics/full_wiot, clear

drop if row_item != 69
drop if col_item > 35

keep year value col_item col_country
rename value GO_usd
rename col_item rnr

joinby rnr using Numerics/SectorsWIOD
keep year col_country GO_usd rnr code sector
rename col_country iso

*weights within country but for different rnr contributing to a sector
bysort iso year sector: egen GO_tot = total(GO_usd)
gen w_raw = GO_usd/GO_tot

* Now renormalize weights within iso x sector4 to sum to 1 (robust to complete missing)
bysort iso year sector: egen wsum = total(w_raw)
gen w = w_raw/wsum

keep year iso GO_usd rnr code sector w
save Numerics/go_sectoral_weights, replace

use Numerics/go_sectoral_weights, clear
drop w
*weights of each country and its each sector into world go
bysort year: egen GO_tot = total(GO_usd)
gen w_raw = GO_usd/GO_tot

* Now renormalize weights within iso x sector4 to sum to 1 (robust to complete missing)
bysort year: egen wsum = total(w_raw)
gen w = w_raw/wsum

drop w_raw
drop wsum
save Numerics/go_country_weights_raw, replace

collapse (sum) GO_usd w, by(year iso sector)
bysort year: egen wsum = total(w)
gen w_n = w/wsum
drop w_n
drop wsum
save Numerics/go_country_weights_sectoral, replace

collapse (sum) GO_usd w, by(year iso)
save Numerics/go_country_weights, replace


*************** INPUT EXCHANGE RATES ******************************************

import delimited Numerics/Input_Files/exchange.csv, clear
rename country country_i
replace country_i = lower(country_i)
rename iso iso_i
rename er er_i
save Numerics/er_i, replace

use Numerics/er_i, clear
rename iso_i iso
joinby year iso using Numerics/go_country_weights
generate row_region = "RoW"
replace row_region = "USA" if iso == "USA"
replace row_region = "CHN" if iso == "CHN"
replace row_region = "IND" if iso == "IND"
replace row_region = "MEX" if iso == "MEX"
replace row_region = "JPN" if iso == "JPN"
replace row_region = "BRA" if iso == "BRA"
replace row_region = "DEU" if iso == "DEU"
replace row_region = "AUS" if iso == "AUS"

drop iso
rename row_region iso
collapse (mean) er_i [aw = GO_usd], by(year iso)
rename iso iso_i
sort iso year

save Numerics/er_short, replace


********** (TAKES LONG TIME TO RUN) Cleaned WIOT 2014 STATA full .dta **********

use Numerics/Input_Files/WIOTS_in_STATA/wiot_full, clear
save Numerics/full_wiot, replace

use Numerics/full_wiot, clear

rename row_item rnr
joinby rnr using Numerics/SectorsWIOD
drop industry
drop code
rename rnr row_item
rename col_item rnr
rename sector row_sector
joinby rnr using Numerics/SectorsWIOD
drop industry
drop code
rename rnr col_item
rename sector col_sector

save Numerics/wiot, replace


**********************BETA_i CALCULATION****************************************

*compensation of employees (in millions of LCU)
import delimited Numerics/Input_Files/comp_i, clear
bysort iso_i code: ipolate comp_i year, gen(interpolated) epolate
bysort iso_i code: replace comp_i = interpolated if year == 2010
bysort iso_i code: replace comp_i = interpolated if year == 2011
drop interpolated
drop if comp_i == .
save Numerics/comp_lcu, replace

joinby code using Numerics/SectorsWIOD
keep iso_i year comp_i sector

generate iso = "RoW"
replace iso = "USA" if iso_i == "USA"
replace iso = "CHN" if iso_i == "CHN"
replace iso = "IND" if iso_i == "IND"
replace iso = "MEX" if iso_i == "MEX"
replace iso = "JPN" if iso_i == "JPN"
replace iso = "BRA" if iso_i == "BRA"
replace iso = "DEU" if iso_i == "DEU"
replace iso = "AUS" if iso_i == "AUS"
drop iso_i
rename iso iso_i

collapse (sum) comp_i, by(iso_i year sector)
save Numerics/comp_lcu_short, replace

reshape wide comp_i, i(year iso_i) j(sector) string
rename comp_iG compG_i
rename comp_iH compH_i
rename comp_iL compL_i
save Numerics/comp_lcu_sectoral, replace

joinby iso_i year using Numerics/er_short
gen compG_i_usd = compG_i * er_i
gen compH_i_usd = compH_i * er_i
gen compL_i_usd = compL_i * er_i

save Numerics/comp_sectoral_usd, replace

**** INT ii_G_i ii_i, VA va_G_i, GO Y_G_i DATA from WIOT (in mn of US $) *******

use Numerics/wiot, clear

drop if row_item < 36
drop if missing(row_sector)

*ADDON: **Shorter Dataset***
generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country
rename col_region col_country

/*
i is importing (so expense of i on intermediates) - i uses INT then adds VA to generate GO which is used as intermediate inputs, final consumption or investment
*/

drop if inlist(col_sector, "CON", "INV")
collapse (sum) value, by(row_sector year col_country col_sector)
rename col_country iso_i
rename row_sector variable
rename col_sector sector
reshape wide value, i(iso_i year sector) j(variable) string
rename valueGO Y_sectoral_i
rename valueINT ii_sectoral_i
rename valueVA va_sectoral_i

save Numerics/wiot_ii_va_go, replace

reshape wide Y_sectoral_i ii_sectoral_i va_sectoral_i, i(iso_i year) j(sector) string
rename Y_sectoral_iG Y_G_usd_i
rename Y_sectoral_iH Y_H_usd_i
rename Y_sectoral_iL Y_L_usd_i
rename ii_sectoral_iG ii_G_usd_i
rename ii_sectoral_iH ii_H_usd_i
rename ii_sectoral_iL ii_L_usd_i
rename va_sectoral_iG va_G_usd_i
rename va_sectoral_iH va_H_usd_i
rename va_sectoral_iL va_L_usd_i

save Numerics/ii_va_go_usd_sectoral, replace
save Numerics/wiot_ii_va_go, replace


*alpha (marla) = betaK

use Numerics/comp_sectoral_usd, clear
joinby iso_i year using Numerics/wiot_ii_va_go

gen betaKG_i = (va_G_usd_i-compG_i_usd)/va_G_usd_i
gen betaKH_i = (va_H_usd_i-compH_i_usd)/va_H_usd_i
gen betaKL_i = (va_L_usd_i-compL_i_usd)/va_L_usd_i

keep iso_i year betaKG_i betaKH_i betaKL_i
sort iso_i year

save Numerics/betaK_alpha, replace

drop if iso_i == "RoW"
collapse (mean) betaKG_i betaKH_i betaKL_i, by(year)
gen iso_i = "RoW"

save Numerics/betaK_row, replace
use Numerics/betaK_alpha, clear
drop if iso_i == "RoW"
append using Numerics/betaK_row
replace iso_i = "RoW" if iso_i ==""

sort iso_i year
save Numerics/betaK_alpha, replace

*****IO tables IO_GH (inputs used from H to G) and IO_G (total inputs used by G)
**** col_country is importer, row is exporter
 
use Numerics/wiot, clear
drop if row_item > 35
drop if col_item > 35

*ADDON: **Shorter Dataset***
generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country
rename col_region col_country

collapse (sum) value, by(year col_country col_sector row_sector)
***IO_GH means inputs from H to G
reshape wide value, i(year col_country row_sector) j(col_sector) string
rename valueG IO_G
rename valueH IO_H
rename valueL IO_L
reshape wide IO_G IO_H IO_L, i(year col_country) j(row_sector) string

rename col_country iso_i
sort iso_i year

save Numerics/IO_MX_sectoral, replace


********Prices caliberation csv

**GGDC 2005 benchmark

import delimited using GGDC_2005_benchmark_35ind.csv, clear
save Numerics/ggdc_2005, replace

**GO shares from WIOD as weights

use Numerics/full_wiot, clear

keep if year==2005
drop if row_item != 69
drop if col_item > 35

keep value col_item col_country
rename value GO_usd
rename col_item rnr

joinby rnr using Numerics/SectorsWIOD
keep col_country GO_usd code sector
rename col_country iso

joinby iso code using Numerics/ggdc_2005
rename value ggdc_p

save Numerics/raw_ggdc_p_and_go_usd, replace

*weights
bysort iso sector: egen GO_tot = total(GO_usd)
gen w_raw = GO_usd/GO_tot

* Now renormalize weights within iso x sector4 to sum to 1 (robust to complete missing)
*bysort iso sector: egen wsum = total(w_raw)
*gen w = w_raw/wsum

keep iso GO_usd sector w_raw
save Numerics/w_raw_2005, replace

use Numerics/raw_ggdc_p_and_go_usd, clear
bysort iso sector: egen GO_tot = total(GO_usd)
gen w_raw = GO_usd/GO_tot

*since my w sums to 1 so rename w_raw to w_raw
rename w_raw w
gen contrib = w * ggdc_p

* Collapse to one row per iso x sector
collapse (sum) P_2005_ = contrib, by(iso sector)

reshape wide P_2005_, i(iso) j(sector) string

sort iso
save Numerics/ggdc_2005_sectoral, replace

use Numerics/go_country_weights, clear
keep if year == 2005

save Numerics/go_country_weight_2005, replace

generate row_region = "RoW"
replace row_region = "USA" if iso == "USA"
replace row_region = "CHN" if iso == "CHN"
replace row_region = "IND" if iso == "IND"
replace row_region = "MEX" if iso == "MEX"
replace row_region = "JPN" if iso == "JPN"
replace row_region = "BRA" if iso == "BRA"
replace row_region = "DEU" if iso == "DEU"
replace row_region = "AUS" if iso == "AUS"
drop iso
rename row_region iso
collapse (sum) GO_usd, by(year iso)

save Numerics/go_short_weights_2005, replace

use Numerics/ggdc_2005_sectoral, clear

joinby iso using Numerics/go_country_weight_2005

generate row_region = "RoW"
replace row_region = "USA" if iso == "USA"
replace row_region = "CHN" if iso == "CHN"
replace row_region = "IND" if iso == "IND"
replace row_region = "MEX" if iso == "MEX"
replace row_region = "JPN" if iso == "JPN"
replace row_region = "BRA" if iso == "BRA"
replace row_region = "DEU" if iso == "DEU"
replace row_region = "AUS" if iso == "AUS"
drop iso
rename row_region iso
collapse (mean) P_2005_G P_2005_H P_2005_L [aw = GO_usd], by(iso)

save Numerics/ggdc_2005_final, replace

** Changes in prices from nominal and real GO from SEA/ go_ppi

import delimited using Numerics/Input_Files/go_p.csv, clear
rename iso_i iso
joinby iso year code using Numerics/go_country_weights_raw 

destring ppi_i, replace force
sort iso code year
by iso code: ipolate ppi_i year, gen(ppi_interpolated) epolate
by iso code: replace ppi_i = ppi_interpolated if year == 2010
by iso code: replace ppi_i = ppi_interpolated if year == 2011
drop ppi_interpolated

collapse (mean) ppi_i [aw = GO_usd], by(year iso sector)
sort iso sector year


* make short dataset with row
joinby year iso using Numerics/go_country_weights
generate row_region = "RoW"
replace row_region = "USA" if iso == "USA"
replace row_region = "CHN" if iso == "CHN"
replace row_region = "IND" if iso == "IND"
replace row_region = "MEX" if iso == "MEX"
replace row_region = "JPN" if iso == "JPN"
replace row_region = "BRA" if iso == "BRA"
replace row_region = "DEU" if iso == "DEU"
replace row_region = "AUS" if iso == "AUS"

drop iso
rename row_region iso
collapse (mean) ppi_i [aw = GO_usd], by(year iso sector)
sort iso year sector

rename ppi_i p_sea_

reshape wide p_sea_, i(iso year) j(sector) string

save Numerics/p_sea, replace

keep if year == 2005
rename p_sea_G p_sea05_G
rename p_sea_L p_sea05_L
rename p_sea_H p_sea05_H

save Numerics/p_sea_2005, replace

use Numerics/p_sea, clear
joinby iso using Numerics/p_sea_2005

replace p_sea_G = p_sea_G/p_sea05_G
replace p_sea_L = p_sea_L/p_sea05_L
replace p_sea_H = p_sea_H/p_sea05_H

keep year iso p_sea_G p_sea_H p_sea_L
save Numerics/p_sea_final, replace

*** adjusted price = p_sea_* x p_2005_*

use Numerics/p_sea_final, clear
joinby iso using Numerics/ggdc_2005_final

gen p_G = p_sea_G * P_2005_G
gen p_L = p_sea_L * P_2005_L
gen p_H = p_sea_H * P_2005_H

keep iso year p_G p_L p_H

save Numerics/final_price_data, replace


*********sectoral labor (per thousand)

import delimited Numerics/Input_Files/emp_sectoral.csv, clear

reshape wide emp, i(iso year) j(sector) string

save Numerics/emp_sectoral, replace


**** Real capital EPOLATE/IPOLATE 2010 2011

import delimited Numerics/Input_Files/real_cap.csv, clear
*rename v1 iso_i
*rename v2 sector
*rename v3 year 
*rename v4 real_k 
replace real_k = "." if real_k == "NA"

destring real_k, replace force

sort iso_i code year
by iso_i code: ipolate real_k year, gen(k_interpolated) epolate
by iso_i code: replace real_k = k_interpolated if year == 2008
by iso_i code: replace real_k = k_interpolated if year == 2009
by iso_i code: replace real_k = k_interpolated if year == 2010
by iso_i code: replace real_k = k_interpolated if year == 2011
drop k_interpolated

sort iso_i year

joinby code using Numerics/SectorsWIOD
drop industry
drop code
drop rnr

generate iso = "RoW"
replace iso = "USA" if iso_i == "USA"
replace iso = "CHN" if iso_i == "CHN"
replace iso = "IND" if iso_i == "IND"
replace iso = "MEX" if iso_i == "MEX"
replace iso = "JPN" if iso_i == "JPN"
replace iso = "BRA" if iso_i == "BRA"
replace iso = "DEU" if iso_i == "DEU"
replace iso = "AUS" if iso_i == "AUS"
drop iso_i

collapse (sum) real_k, by(iso year sector)
rename real_k real_k_
reshape wide real_k_, i(iso year) j(sector) string
save Numerics/real_k_sectoral_new, replace

**** K_USD = real_K * GFCF_P (k_ppi)

import delimited Numerics/Input_Files/cap_ppi.csv, clear
*rename v1 iso_i
*rename v2 sector
*rename v3 year 
*rename v4 real_k 
replace k_ppi = "." if k_ppi == "NA"

destring k_ppi, replace force

sort iso_i code year
by iso_i code: ipolate k_ppi year, gen(k_interpolated) epolate
by iso_i code: replace k_ppi = k_interpolated if year == 2008
by iso_i code: replace k_ppi = k_interpolated if year == 2009
by iso_i code: replace k_ppi = k_interpolated if year == 2010
by iso_i code: replace k_ppi = k_interpolated if year == 2011
drop k_interpolated

sort iso_i year

joinby code using Numerics/SectorsWIOD
drop industry
drop code
drop rnr

generate iso = "RoW" 
replace iso = "USA" if iso_i == "USA"
replace iso = "CHN" if iso_i == "CHN"
replace iso = "IND" if iso_i == "IND"
replace iso = "MEX" if iso_i == "MEX"
replace iso = "JPN" if iso_i == "JPN"
replace iso = "BRA" if iso_i == "BRA"
replace iso = "DEU" if iso_i == "DEU"
replace iso = "AUS" if iso_i == "AUS"
drop iso_i

joinby year iso sector using Numerics/go_country_weights_sectoral

collapse (mean) k_ppi [aw = GO_usd], by(year iso sector)
rename k_ppi k_ppi_
reshape wide k_ppi_, i(iso year) j(sector) string

joinby iso year using Numerics/real_k_sectoral_new
rename iso iso_i
joinby iso_i year using Numerics/er_short

gen k_G_usd = real_k_G * (k_ppi_G / 100) * er_i
gen k_H_usd = real_k_H * (k_ppi_H / 100) * er_i
gen k_L_usd = real_k_L * (k_ppi_L / 100) * er_i

save Numerics/k_usd_new, replace
 

 ******PC, PX, inv and FD
 
use Numerics/wiot, clear
drop if col_item < 37
drop if row_item > 35

*ADDON: **Shorter Dataset***
generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country
rename col_region col_country

replace col_sector = "PX" if col_item == 41

collapse (sum) value, by(year col_country col_sector row_sector)
rename col_country iso_j
reshape wide value, i(iso_j year row_sector) j(col_sector) string
rename valueCON PC_
rename valueINV INV_
rename valuePX PX_
reshape wide PC_ INV_ PX_, i(iso_j year) j(row_sector) string

save Numerics/PC_PX_INV_usd, replace


*****Bilateral_trade

use Numerics/wiot, clear

drop if row_item > 35
*ADDON: **Shorter Dataset***
generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country
rename col_region col_country

collapse (sum) value, by(year row_country col_country row_sector)
* Each obs is X_{i->j}^s with i=row_country, j=col_country, s=row_sector
*reshape wide value, i(year row_country col_country) j(row_sector) string

gen exporter_num = 1
replace exporter_num = 2 if row_country == "BRA"
replace exporter_num = 3 if row_country == "CHN"
replace exporter_num = 4 if row_country == "DEU"
replace exporter_num = 5 if row_country == "IND"
replace exporter_num = 6 if row_country == "JPN"
replace exporter_num = 7 if row_country == "MEX"
replace exporter_num = 8 if row_country == "RoW"
replace exporter_num = 9 if row_country == "USA"

gen importer_num = 1
replace importer_num = 2 if col_country == "BRA"
replace importer_num = 3 if col_country == "CHN"
replace importer_num = 4 if col_country == "DEU"
replace importer_num = 5 if col_country == "IND"
replace importer_num = 6 if col_country == "JPN"
replace importer_num = 7 if col_country == "MEX"
replace importer_num = 8 if col_country == "RoW"
replace importer_num = 9 if col_country == "USA"

keep year row_sector exporter_num importer_num value
gen sector = .
replace sector = 1 if row_sector == "G"
replace sector = 2 if row_sector == "L"
replace sector = 3 if row_sector == "H"

rename value bt_all
sort importer_num exporter_num year sector

save Numerics/bt_all, replace


******************** IMPORT-EXPORT DATA X_ and M_ ******************************

***exports
use Numerics/wiot, clear
drop if row_item > 35

generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country 
rename col_region col_country 

*row_country exports for both intermediate and final consumption 
collapse (sum) value, by(year row_country row_sector)
reshape wide value, i(year row_country) j(row_sector) string
rename valueG X_i_G
rename valueH X_i_H
rename valueL X_i_L
rename row_country iso_i

sort iso_i year

save Numerics/X_i, replace


***imports
use Numerics/wiot, clear
drop if row_item > 35

generate row_region = "RoW"
replace row_region = "USA" if row_country == "USA"
replace row_region = "CHN" if row_country == "CHN"
replace row_region = "IND" if row_country == "IND"
replace row_region = "MEX" if row_country == "MEX"
replace row_region = "JPN" if row_country == "JPN"
replace row_region = "BRA" if row_country == "BRA"
replace row_region = "DEU" if row_country == "DEU"
replace row_region = "AUS" if row_country == "AUS"

generate col_region = "RoW"
replace col_region = "USA" if col_country == "USA"
replace col_region = "CHN" if col_country == "CHN"
replace col_region = "IND" if col_country == "IND"
replace col_region = "MEX" if col_country == "MEX"
replace col_region = "JPN" if col_country == "JPN"
replace col_region = "BRA" if col_country == "BRA"
replace col_region = "DEU" if col_country == "DEU"
replace col_region = "AUS" if col_country == "AUS"

drop row_country
drop col_country
rename row_region row_country 
rename col_region col_country 

*col_country imports row_sector for both intermediate and final consumption 
collapse (sum) value, by(year col_country row_sector)
reshape wide value, i(year col_country) j(row_sector) string
rename valueG M_i_G
rename valueH M_i_H
rename valueL M_i_L
rename col_country iso_i

sort iso_i year

save Numerics/M_i, replace





*export delimited using "sector_price_paths_rel_2000_2014.csv", replace



