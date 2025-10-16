clear all
do Globals

import delimited "${Data}/StateAnalysisFileTfp.csv", clear case(preserve)

loc samp inrange(year, 1994, 2021)
loc sampno21 inrange(year, 1994, 2020)

loc graphs "Wage00Response_Iv1990 Wage01Response_Iv1990 Wage00Response_Iv1990_F Wage01Response_Iv1990_F"

/*************************************************
Some Cleaning - 
Identical to the cleaning part in the Z regs code
*************************************************/

* Change statefip to string of length 2
tostring statefip, replace
replace statefip = "0" + statefip if strlen(statefip) == 1

* Total employment
gen emp = BodiesSupplied00 + BodiesSupplied01

* Generate foreign-born labor share and its growth rate
gen f = BodiesSupplied01
bys statefip (year): gen fg = (f - f[_n - 1]) / emp

/* The current dataset does not have the labor supplies of each 
migrant group over time. It only has 1990 and prior. So, let's
fetch these.
*/
frame create GroupBreakdown
frame GroupBreakdown {

    * Recall that, before we estimated TFP, this version of the file was in long format and had information
    * on the labor supplies of each migrant group
    use "${Data}/StateAnalysisFile.dta", clear
    replace foreign = 1 if ImmigrantGroup != "United States"
    keep ImmigrantGroup statefip year BodiesSupplied StateName

    * Prepare to reshape wide
    gen abbr = "_CaAuNz" if ImmigrantGroup == "Canada-Australia-New Zealand"
    replace abbr = "_US" if ImmigrantGroup == "United States"
    replace abbr = "_WestEu" if ImmigrantGroup == "Western Europe"
    replace abbr = "_LA" if ImmigrantGroup == "Latin America"
    replace abbr = "_EastEu" if ImmigrantGroup == "Russia and Eastern Europe"
    replace abbr = "_AsiaOther" if ImmigrantGroup == "Rest of Asia"
    replace abbr = "_China" if ImmigrantGroup == "China"
    replace abbr = "_Africa" if ImmigrantGroup == "Africa"
    replace abbr = "_Mexico" if ImmigrantGroup == "Mexico"
    replace abbr = "_Other" if ImmigrantGroup == "Other"
    replace abbr = "_India" if ImmigrantGroup == "India"

    drop ImmigrantGroup StateName
    
    * Will have to come back to figure out what is happening with the duplicates
    * For now, keep the duplicates which is a larger part of the foreign born
    gsort abbr statefip year -BodiesSupplied
    duplicates drop abbr statefip year, force

    reshape wide Bodies*, i(statefip year) j(abbr) s

    qui ds Bodies*
    foreach v in `r(varlist)' {
        replace `v' = 0 if mi(`v')
    }

}

frlink m:1 statefip year, frame(GroupBreakdown)
frget Bodies* , from(GroupBreakdown)
drop GroupBreakdown

* Calculate the shares
qui ds BodiesSupplied_*
foreach v in `r(varlist)' {

    loc region = subinstr("`v'", "BodiesSupplied_", "", 1)

    * Share of the foreign born population
    if "`region'" != "US" {

        gen s_`region' = `v' / emp
        bys statefip (year): gen s_`region'_L1 = s_`region'[_n - 1]
        bys statefip (year): gen s_`region'_L2 = s_`region'[_n - 2]
    }

}

* Create fixed 1990 shares
qui ds *1990
loc vars1990 "`r(varlist)'"
egen emp1990 = rowtotal(`vars1990')
gen BodiesSupplied01_1990 = emp1990 - US1990

foreach v in `vars1990' {
    
    loc region = subinstr("`v'", "1990", "", 1)

    if "`region'" != "US" gen s_`region'_1990 = `v' / emp1990

}

* Create aggregate shifts. Can use the GroupBreakdown frame again
frame GroupBreakdown {

    qui ds statefip year, not
    collapse (sum) `r(varlist)', by(year)
    egen emp = rowtotal(Bodies*)
    qui ds Bodies*
    foreach v in `r(varlist)' {
        loc region = subinstr("`v'", "BodiesSupplied_", "", 1)
        gen fg_agg_`region' = (`v' - `v'[_n - 1]) / emp 
    }

    replace fg_agg_Africa = . if year == 1996 // outlier
}

frlink m:1 year, frame(GroupBreakdown)
frget fg_agg*, from(GroupBreakdown)
drop GroupBreakdown


* Create the Bartik instruments
qui ds fg_agg_*
foreach v in `r(varlist)' {

    loc region = subinstr("`v'", "fg_agg_", "", 1)
    if !inlist("`region'", "Other", "US") {
        gen Bartik_1990_`region' = s_`region'_1990 * `v'
        gen Bartik_L1_`region' = s_`region'_L1 * `v'
        gen Bartik_L2_`region' = s_`region'_L2 * `v'
    }

}

egen Bartik_1990 = rowtotal(Bartik_1990_*), missing // Pre-period shares
egen Bartik_L1 = rowtotal(Bartik_L1_*), missing     // Lagged shares
egen Bartik_L2 = rowtotal(Bartik_L2_*), missing

* Calculate task aggregate growth rates
forvalues h = -3/9 { // LHS variables
    if `h' < 0 loc name = "L" + string(abs(`h'))
    if `h' >= 0 loc name = "F" + string(abs(`h'))
    bys statefip (year): gen W01g_`name' = Wage01[_n + `h']/Wage01[_n - 1] - 1 // NOTE: These were already deflated when I estimated TFP
    bys statefip (year): gen W00g_`name' = Wage00[_n + `h']/Wage00[_n - 1] - 1
}

encode statefip, gen(state)
xtset state year

/*************************************************
REGRESSIONS
*************************************************/

* Fixed-shares design
frame create Estimates
frame Estimates {

    gen h = .
    gen BetaForeignIv1990 = .
    gen SeForeignIv1990 = .
    gen FForeignIv1990 = .

    gen BetaDomesticIv1990 = .
    gen SeDomesticIv1990 = .
    gen FDomesticIv1990 = .

}

forvalues h = -3/9 {

    if `h' != -1 {
        di "***********************************************************************************************************"
        di "COMPUTING REGRESSION AT HORIZON `h'"
        di "***********************************************************************************************************"
        if `h' < 0 loc horizon = "L" + string(abs(`h'))
        if `h' >= 0 loc horizon = "F" + string(abs(`h'))

        * Foreign born
        qui ivreghdfe W01g_`horizon' (fg = Bartik_1990) [pw = emp] if `samp', absorb(state year) vce(robust)

        * Record the results in the Estimates frame
        frame Estimates {
            insobs 1
            replace h = `h' if _n == _N
            replace BetaForeignIv1990 = _b[fg] if _n == _N
            replace SeForeignIv1990 = _se[fg] if _n == _N
            replace FForeignIv1990 = `e(widstat)' if _n == _N
        }

        * Domestic
        qui ivreghdfe W00g_`horizon' (fg = Bartik_1990) [pw = emp] if `samp', absorb(state year) vce(robust)

        frame Estimates {
            insobs 1
            replace h = `h' if _n == _N
            replace BetaDomesticIv1990 = _b[fg] if _n == _N
            replace SeDomesticIv1990 = _se[fg] if _n == _N
            replace FDomesticIv1990 = `e(widstat)' if _n == _N
        }
    
    }

    * Fill in the horizon -1 results
    else {
        frame Estimates {
            insobs 1
            replace h = `h' if _n == _N
            
            replace BetaForeignIv1990 = 0 if _n == _N
            replace SeForeignIv1990 = . if _n == _N
            replace FForeignIv1990 = . if _n == _N

            replace BetaDomesticIv1990 = 0 if _n == _N
            replace SeDomesticIv1990 = . if _n == _N
            replace FDomesticIv1990 = . if _n == _N
        }
    }
    
}

/*************************************************
GRAPHS
*************************************************/

* IRF of pre-period share
frame Estimates {

    * Confidence intervals - 95%
    gen BetaForeignIv1990_upper = BetaForeignIv1990 + 1.96 * SeForeignIv1990
    gen BetaForeignIv1990_lower = BetaForeignIv1990 - 1.96 * SeForeignIv1990

    gen BetaDomesticIv1990_upper = BetaDomesticIv1990 + 1.96 * SeDomesticIv1990
    gen BetaDomesticIv1990_lower = BetaDomesticIv1990 - 1.96 * SeDomesticIv1990

    tw connected BetaForeignIv1990 h if inrange(h,-3, 9), ms(oh) mc("0 147 245") xlab(-3(1)9, nogrid) sort || rcap BetaForeignIv1990_upper BetaForeignIv1990_lower h, lcolor("0 147 245") ylab(, nogrid) ///
    ytitle("{&epsilon}{subscript:F}") xtitle("Horizon") legend(off) yline(0, lc(black%50) lp(solid)) name(Wage01Response_Iv1990)

    tw connected BetaDomesticIv1990 h if inrange(h,-3, 9), ms(oh) mc("0 147 245") xlab(-3(1)9, nogrid) sort || rcap BetaDomesticIv1990_upper BetaDomesticIv1990_lower h, lcolor("0 147 245") ylab(, nogrid) ///
    ytitle("{&epsilon}{subscript:D}") xtitle("Horizon") legend(off) yline(0, lc(black%50) lp(solid)) name(Wage00Response_Iv1990)

}

* F Stats of IRF for pre-period
frame Estimates {

    graph bar FForeignIv1990 if h > -1, over(h) bar(1, color("0 147 245") fcolor("0 147 245")) ylab(0(20)180, nogrid labsize(small)) ///
    yline(10,lc(black%70) lp(dash)) legend(off) b1title("Horizon") ytitle("First Stage F") name(Wage01Response_Iv1990_F)

    graph bar FDomesticIv1990 if h > -1, over(h) bar(1, color("0 147 245") fcolor("0 147 245")) ylab(0(20)180, nogrid labsize(small)) ///
    yline(10,lc(black%70) lp(dash)) legend(off) b1title("Horizon") ytitle("First Stage F") name(Wage00Response_Iv1990_F)
    
}

* Export
foreach g in `graphs' {

    graph export "${Graphs}/`g'.pdf", name(`g') replace

}
