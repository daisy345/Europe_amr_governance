* HAMR staggered-DID event study: full observable event-time window.
* Third palette trial: named colors red / orange / navy / emerald (CS / Gardner / Stacked / TWFE).

capture log close _all
clear all
version 18.0
set more off
local projectdir "`c(pwd)'"
adopath ++ "`projectdir'/code/vendor/csdid"
adopath ++ "`projectdir'/code/vendor/drdid"
local datadir "`projectdir'/input"
local outdir "`projectdir'/output"
log using "`outdir'/Human_AMR_score.log", replace text


foreach cmd in ftools reghdfe did2s event_plot csdid stackedev {
    cap which `cmd'
    if _rc {
        di as error "Required command `cmd' is unavailable."
        exit 199
    }
}

use "`datadir'/Human_AMR_score.dta", clear
tsset i t

* The regenerated dataset ships Ei, K, gvar, treat_year, no_treat and the
* event-time dummies F*event/L*event prebuilt (K=-1 omitted); reuse them.
assert treat_year == Ei

* Determine the complete observable event-time support from the data.
quietly summarize K, meanonly
local kmin = r(min)
local kmax = r(max)
di as result "Observable event-time range: `kmin' to `kmax'"

* 1. Callaway and Sant'Anna (2020)
csdid Y, ivar(i) time(t) gvar(gvar) agg(event)
matrix cs_b = e(b)
matrix cs_v = e(V)

* 2. Gardner DID-2S (2021)
did2s Y, first_stage(i.i i.t) second_stage(F*event L*event) treatment(D) cluster(i)
matrix did2s_b = e(b)
matrix did2s_v = e(V)

* 3. Stacked event-study DID (Cengiz et al., 2019)
foreach var of varlist F*event L*event {
    replace `var' = 0 if no_treat
}
preserve
stackedev Y F*event L*event, cohort(treat_year) time(t) ///
    never_treat(no_treat) unit_fe(i) clust_unit(i)
restore
matrix stackedev_b = e(b)
matrix stackedev_v = e(V)

* 4. Conventional TWFE event study (reference: K=-1)
reghdfe Y F*event L*event, absorb(i t) vce(cluster i)
estimates store ols

* Combined figure: third palette trial (red/orange/navy/emerald), vsmall
* filled-circle markers, thin 60%-opacity CI lines with short caps,
* series dodged by -0.25(0.09)0.25.
event_plot cs_b#cs_v did2s_b#did2s_v stackedev_b#stackedev_v ols, ///
    stub_lag(Tp# L#event L#event L#event) ///
    stub_lead(Tm# F#event F#event F#event) ///
    plottype(scatter) ciplottype(rcap) together ///
    perturb(-0.25(0.09)0.25) noautolegend ///
    trimlead(`=abs(`kmin')') trimlag(`kmax') ///
    graph_opt( ///
    xtitle("Time since NAP launch (years)") ///
    ytitle("Human AMR score") xlabel(`kmin'(1)`kmax', nogrid) ///
    xscale(range(`kmin' `kmax') extend) ylabel(, angle(horizontal) nogrid) ///
    legend(order(1 "Callaway-Sant'Anna" 3 "Gardner DID-2S" 5 "Stacked DID" 7 "TWFE OLS") ///
    rows(2) position(6) region(style(none))) ///
    xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) ///
    graphregion(color(white)) bgcolor(white)) ///
    /// 1) Callaway-Sant'Anna in red
    lag_opt1(mcolor(red) lcolor(red) msymbol(circle) msize(vsmall)) ///
    lead_opt1(mcolor(red) lcolor(red) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt1(lcolor(red%60) lwidth(thin) msize(small)) ///
    lead_ci_opt1(lcolor(red%60) lwidth(thin) msize(small)) ///
    /// 2) Gardner DID-2S in orange
    lag_opt2(mcolor(orange) lcolor(orange) msymbol(circle) msize(vsmall)) ///
    lead_opt2(mcolor(orange) lcolor(orange) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt2(lcolor(orange%60) lwidth(thin) msize(small)) ///
    lead_ci_opt2(lcolor(orange%60) lwidth(thin) msize(small)) ///
    /// 3) Stacked DID in navy
    lag_opt3(mcolor(navy) lcolor(navy) msymbol(circle) msize(vsmall)) ///
    lead_opt3(mcolor(navy) lcolor(navy) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt3(lcolor(navy%60) lwidth(thin) msize(small)) ///
    lead_ci_opt3(lcolor(navy%60) lwidth(thin) msize(small)) ///
    /// 4) TWFE OLS in emerald
    lag_opt4(mcolor(emerald) lcolor(emerald) msymbol(circle) msize(vsmall)) ///
    lead_opt4(mcolor(emerald) lcolor(emerald) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt4(lcolor(emerald%60) lwidth(thin) msize(small)) ///
    lead_ci_opt4(lcolor(emerald%60) lwidth(thin) msize(small))
* width(2155): match the resolution of the GUI-exported v2 PNGs.
graph export "`outdir'/Human_AMR_score.png", width(2155) replace
graph export "`outdir'/Human_AMR_score.pdf", replace

log close
