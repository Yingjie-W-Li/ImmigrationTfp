********************************************************************************
* WageRegressions.do
* Local-Projection IV (LPIV) for (foreign-born) wages using Bartik instruments
* - Main IV: fixed pre-period (1994) shares × leave-one-out national shocks
* - Alt  IV: lagged shares × leave-one-out national shocks
* Data files used:
*   - DATA/StateAnalysisFile.dta     (state × immigrant-group × year panel)
*   - DATA/StateAnalysisFileTfp.dta  (state-year aggregates incl. BodiesSupplied0/1)
********************************************************************************

version 17
clear all
set more off
set matsize 11000

* === Paths ===
do Globals

********************************************************************************
* Build state-year outcomes: foreign-born wage (state-year avg)
*    and endogenous migration flow f_t
********************************************************************************

* Foreign-born state-year wage
frame create Outcome
frame Outcome: use "$Data/StateAnalysisFile.dta", clear
frame Outcome: keep if inrange(year, 1989, .)
frame Outcome: keep ImmigrantGroup foreign statefip StateName year Wage BodiesSupplied

frame Outcome: bysort statefip year: egen BodiesForeign_sy = total(cond(foreign==1, BodiesSupplied, .))
frame Outcome: gen w_contrib = cond(foreign==1, Wage * BodiesSupplied, .)
frame Outcome: bysort statefip year: egen w_num_sy = total(w_contrib)

frame Outcome: gen wF = w_num_sy / BodiesForeign_sy
frame Outcome: drop w_contrib w_num_sy

frame Outcome: keep if !missing(wF)
frame Outcome: bysort statefip year: keep if _n==1
frame Outcome: keep statefip StateName year wF
frame Outcome: gen logwF = log(wF)

frame change Outcome
frame put statefip StateName year logwF, into(Outcome_sy)

* Migration flow f_t
frame create Agg
frame Agg: use "$Data/StateAnalysisFileTfp.dta", clear
frame Agg: destring year, replace force
frame Agg: encode statefip, gen(StateID)
frame Agg: xtset StateID year
frame Agg: gen F_lead = F.BodiesSupplied1
frame Agg: gen wt     = BodiesSupplied1 + BodiesSupplied0
frame Agg: gen f      = (F_lead - BodiesSupplied1) / wt
frame Agg: keep statefip year wt f

frame change Outcome_sy
frlink 1:1 statefip year, frame(Agg) gen(_ag)
frget wt f, from(_ag)
drop _ag

drop if missing(logwF) | missing(f) | missing(wt)

********************************************************************************
* Bartik IVs from group-level panel:
*    - Fixed 1994 shares (foreign only)
*    - Lagged shares (foreign only)
*    - Leave-one-out national growth shocks by group-year
********************************************************************************

* Shares frame (state × group × year), foreign only
frame create Shares
frame Shares: use "$Data/StateAnalysisFile.dta", clear
frame Shares: keep ImmigrantGroup foreign statefip year BodiesSupplied
frame Shares: keep if foreign==1
frame Shares: destring year, replace force

* Ensure uniqueness of (state,group,year)
frame Shares: isid statefip ImmigrantGroup year
if _rc {
    frame Shares: collapse (sum) BodiesSupplied, by(statefip ImmigrantGroup year)
}

* Fixed pre-period shares 1994
frame Shares: preserve
    frame Shares: keep if year==1994
    frame Shares: bysort statefip: egen BodiesForeign_1994 = total(BodiesSupplied)
    frame Shares: gen share1994 = BodiesSupplied / BodiesForeign_1994
    frame Shares: keep statefip ImmigrantGroup share1994
    frame Shares: tempfile Share1994
    frame Shares: save `Share1994', replace
frame Shares: restore

* Lagged shares s_{sg,t-1}
frame Shares: bysort statefip year: egen BodiesForeign_sy = total(BodiesSupplied)
frame Shares: egen panel_id = group(statefip ImmigrantGroup)
frame Shares: xtset panel_id year
frame Shares: gen BodiesSupplied_l1   = L.BodiesSupplied
frame Shares: bysort statefip year: egen BodiesForeign_sy_l1 = total(BodiesSupplied_l1)
frame Shares: gen share_lag = BodiesSupplied_l1 / BodiesForeign_sy_l1

*  Shifts (national + leave-one-out growth by group-year)
frame create Shifts
frame Shifts: use "$Data/StateAnalysisFile.dta", clear
frame Shifts: keep ImmigrantGroup foreign statefip year BodiesSupplied
frame Shifts: keep if foreign==1
frame Shifts: destring year, replace force

* Make unique (group,state,year)
frame Shifts: collapse (sum) BodiesSupplied, by(ImmigrantGroup statefip year)

* National totals by (group,year), then lags in a clean (group,year) panel
frame Shifts: bysort ImmigrantGroup year: egen NatBodies_gyr = total(BodiesSupplied)
frame Shifts: preserve
    frame Shifts: keep ImmigrantGroup year NatBodies_gyr
    frame Shifts: bysort ImmigrantGroup year: keep if _n==1
    frame Shifts: egen g = group(ImmigrantGroup)
    frame Shifts: xtset g year
    frame Shifts: gen NatBodies_gyr_l1 = L.NatBodies_gyr
    frame Shifts: gen nat_growth_gyr   = (NatBodies_gyr - NatBodies_gyr_l1) / NatBodies_gyr_l1
    tempfile NATIONAL
    frame Shifts: save `NATIONAL', replace
frame Shifts: restore

* Leave-one-out totals by (group,state,year), then lags in (group×state) panels
frame Shifts: gen LoBodies_gyr = NatBodies_gyr - BodiesSupplied
frame Shifts: egen gs = group(ImmigrantGroup statefip)
frame Shifts: xtset gs year
frame Shifts: gen LoBodies_gyr_l1 = L.LoBodies_gyr
frame Shifts: gen loo_growth_gyr  = (LoBodies_gyr - LoBodies_gyr_l1) / LoBodies_gyr_l1

* Reduce LOO growth to (group,year) and add national growth (optional)
frame Shifts: keep ImmigrantGroup year loo_growth_gyr
frame Shifts: bysort ImmigrantGroup year: keep if _n==1
frame Shifts: merge 1:1 ImmigrantGroup year using `NATIONAL', nogenerate

* Attach shocks to Shares; build Bartik pieces and collapse to state-year IVs
frame change Shares
* link (m:1) by (group,year) from Shifts
frlink m:1 ImmigrantGroup year, frame(Shifts) gen(_sh)
frget loo_growth_gyr nat_growth_gyr, from(_sh)
drop _sh

* add fixed pre-period shares
frame Shares: merge m:1 statefip ImmigrantGroup using `Share1994', nogenerate

frame Shares: gen B1994_piece = share1994 * loo_growth_gyr
frame Shares: gen Blag_piece  = share_lag  * loo_growth_gyr

frame Shares: collapse (sum) Bartik1994 = B1994_piece ///
                              (sum) BartikLag  = Blag_piece, by(statefip year)

********************************************************************************
* Merge IVs into state-year outcome panel
********************************************************************************
frame change Outcome_sy
frlink 1:1 statefip year, frame(Shares) gen(_sh2)
frget Bartik1994 BartikLag, from(_sh2)
drop _sh2

drop if missing(Bartik1994) | missing(BartikLag)

********************************************************************************
* Local-Projection setup: build forward changes of log wages
*    y_{t+h} - y_t for h ∈ {0,…,10}
********************************************************************************
encode StateName, gen(State)
destring year, replace force
xtset State year

tempvar basey
gen `basey' = logwF
local Hlist 0 1 2 3 4 5 6 7 8 9 10

foreach h of local Hlist {
    gen dlogwF_f`h' = F`h'.logwF - `basey'
}

label var wt "pweight: total employment (native+foreign)"
label var f  "immigration flow share (F_{t+1}-F_t)/Total"

********************************************************************************
* Estimation loop over horizons
********************************************************************************
capture frame drop Estimates
frame create Estimates str12 spec byte h double b se fs r2 n

local FE i.State i.year
local W  [pw=wt]
local CL , vce(cluster year)

frame change Outcome_sy
postutil clear
tempname mem
postfile `mem' str12 spec byte h double(b se fs r2 n) using "$Out/Wage_LPIV_Results.dta", replace

* skip h==0 in the loop because dlogwF_f0 is always constant and zero
foreach h of local Hlist {
    if `h'==0 continue

    local lhs dlogwF_f`h'

    qui reghdfe `lhs' f `W', absorb(State year) vce(cluster year)
    post `mem' ("OLS") (`h') (_b[f]) (_se[f]) (.) (e(r2)) (e(N))

    qui ivregress 2sls `lhs' (f = Bartik1994) `FE' `W' `CL'
    qui estat firststage
    scalar _fs = .
    capture matrix FS = r(F)
    if !_rc capture scalar _fs = FS[1,1]
    post `mem' ("IV_FPPS") (`h') (_b[f]) (_se[f]) (_fs) (e(r2)) (e(N))

    qui ivregress 2sls `lhs' (f = BartikLag) `FE' `W' `CL'
    qui estat firststage
    scalar _fs2 = .
    capture matrix FS = r(F)
    if !_rc capture scalar _fs2 = FS[1,1]
    post `mem' ("IV_Lagged") (`h') (_b[f]) (_se[f]) (_fs2) (e(r2)) (e(N))
}

postclose `mem'
frame Estimates: use "$Out/Wage_LPIV_Results.dta", clear

/*
********************************************************************************
* Plot β(h) with 90% CIs: OLS vs IV_FPPS vs IV_Lagged
********************************************************************************
frame change Estimates
preserve
keep if inlist(spec, "OLS", "IV_FPPS", "IV_Lagged")
reshape wide b se fs r2 n, i(h) j(spec) string

gen ci_lo_OLS    = bOLS      - 1.64*seOLS
gen ci_hi_OLS    = bOLS      + 1.64*seOLS
gen ci_lo_IV1994 = bIV_1994  - 1.64*seIV_1994
gen ci_hi_IV1994 = bIV_1994  + 1.64*seIV_1994

twoway
 (rarea ci_lo_OLS ci_hi_OLS h, sort lwidth(none) fintensity(20)) ///
 (line  bOLS h, lpattern(solid) lwidth(medthick)) ///
 (rarea ci_lo_IV1994 ci_hi_IV1994 h, sort lwidth(none) fintensity(20)) ///
 (line  bIV_1994 h, lpattern(dash) lwidth(medthick)), ///
 legend(order(2 "OLS β(h)" 4 "IV (1994 shares) β(h)") pos(6) ring(0)) ///
 ytitle("Effect on log w^F_{t+h} - log w^F_t") xtitle("Horizon h") ///
 title("Local Projections: Foreign-born Wages") name(WageLP, replace)

graph export "$Graphs/WageLP_OLS_vs_IVs.pdf", replace
restore
*/

********************************************************************************
* Save xlsx
********************************************************************************
frame Estimates: export excel using "$Out/Wage_LPIV_Results.xlsx", firstrow(variables) replace

********************************************************************************


********************************************************************************
* build wage-incidence elasticities ε_w^F(h) = β_h * \bar f
*     - \bar f is the *level* migrant share (foreign / total), weighted by wt
*     - For each spec & horizon, compute: eps(h) = b(h) * fbar, se_eps(h) = se(h) * fbar
********************************************************************************

capture noisily log close _all

*-- Compute fbar (weighted mean level migrant share)
preserve
use "$Data/StateAnalysisFileTfp.dta", clear
gen wt = BodiesSupplied1 + BodiesSupplied0
gen f_level = BodiesSupplied1 / wt
quietly summarize f_level [aw=wt] if !missing(f_level) & f_level>=0 & f_level<=1
scalar fbar = r(mean)
display as result "Weighted mean migrant share (fbar) = " %9.6f fbar
restore

*-- Load the saved LP results and compute ε_w^F(h)
use "$Out/Wage_LPIV_Results.dta", clear
gen eps   = b  * fbar
gen se_eps = se * fbar

* Keep common specs
keep if inlist(spec, "OLS", "IV_FPPS", "IV_Lagged")

* Save a long-form elasticities table
order spec h b se eps se_eps fs r2 n
label var eps    "epsilon_wF(h) = beta(h) * fbar"
label var se_eps "SE[epsilon_wF(h)] = fbar * SE[beta(h)]"
save "$Out/Wage_LPIV_Elasticities.dta", replace
export excel using "$Out/Wage_LPIV_Elasticities.xlsx", firstrow(variables) replace
