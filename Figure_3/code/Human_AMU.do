* HAMU staggered-DID event study: CS, DID-2S and TWFE on K=-12 to K=13;
* Stacked DID with tails binned to K=-10 through K=10.
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
log using "`outdir'/Human_AMU.log", replace text
foreach cmd in ftools reghdfe did2s event_plot csdid stackedev {
    cap which `cmd'
    if _rc exit 199
}

use "`datadir'/Human_AMU.dta", clear
rename amu Y
rename id i
tsset i t
by i (t): egen Ei = min(cond(D == 1, t, .))
gen K = t - Ei

* CS/DID-2S/TWFE: raw K=-12 through K=13; K=-1 omitted.
forvalues k = -12/-2 {
    local a = abs(`k')
    gen F`a'event = (K == `k')
}
forvalues k = 0/13 {
    gen L`k'event = (K == `k')
}

gen gvar = cond(missing(Ei), 0, Ei)
csdid Y, ivar(i) time(t) gvar(gvar) agg(event)
matrix cs_b = e(b)
matrix cs_v = e(V)

did2s Y, first_stage(i.i i.t) second_stage(F*event L*event) treatment(D) cluster(i)
matrix did2s_b = e(b)
matrix did2s_v = e(V)

* Stacked DID: bin both tails at -10 and +10; K=-1 omitted.
gen K_stacked = K
replace K_stacked = -10 if K < -10
replace K_stacked = 10 if K > 10
forvalues k = -10/-2 {
    local a = abs(`k')
    gen SF`a'event = (K_stacked == `k')
}
forvalues k = 0/10 {
    gen SL`k'event = (K_stacked == `k')
}
gen treat_year = Ei
gen no_treat = missing(Ei)
foreach var of varlist SF*event SL*event {
    replace `var' = 0 if no_treat
}
preserve
stackedev Y SF*event SL*event, cohort(treat_year) time(t) ///
    never_treat(no_treat) unit_fe(i) clust_unit(i)
restore
matrix stacked_b = e(b)
matrix stacked_v = e(V)

reghdfe Y F*event L*event, absorb(i t) vce(cluster i)
estimates store ols

* Combined figure: third palette trial (red/orange/navy/emerald), vsmall
* filled-circle markers, thin 60%-opacity CI lines with short caps,
* series dodged by -0.25(0.08)0.25.
event_plot cs_b#cs_v did2s_b#did2s_v stacked_b#stacked_v ols, ///
    stub_lag(Tp# L#event SL#event L#event) ///
    stub_lead(Tm# F#event SF#event F#event) ///
    plottype(scatter) ciplottype(rcap) together ///
    perturb(-0.25(0.08)0.25) noautolegend trimlead(12) trimlag(13) ///
    graph_opt( ///
    xtitle("Time since NAP launch (years)") ///
    ytitle("Human AMU") xlabel(-12(1)13, nogrid) ///
    xscale(range(-12 13) extend) ylabel(, angle(horizontal) nogrid) ///
    legend(order(1 "Callaway-Sant'Anna" 3 "Gardner DID-2S" 5 "Stacked DID" 7 "TWFE OLS") ///
    rows(2) position(6) region(style(none))) ///
    xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) ///
    graphregion(color(white)) bgcolor(white)) ///
    /// 1) Callaway-Sant'Anna in red
    lag_opt1(mcolor(red) lcolor(red) msymbol(circle) msize(vsmall)) ///
    lead_opt1(mcolor(red) lcolor(red) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt1(lcolor(red%60) lwidth(thin) msize(vsmall)) ///
    lead_ci_opt1(lcolor(red%60) lwidth(thin) msize(vsmall)) ///
    /// 2) Gardner DID-2S in orange
    lag_opt2(mcolor(orange) lcolor(orange) msymbol(circle) msize(vsmall)) ///
    lead_opt2(mcolor(orange) lcolor(orange) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt2(lcolor(orange%60) lwidth(thin) msize(vsmall)) ///
    lead_ci_opt2(lcolor(orange%60) lwidth(thin) msize(vsmall)) ///
    /// 3) Stacked DID in navy
    lag_opt3(mcolor(navy) lcolor(navy) msymbol(circle) msize(vsmall)) ///
    lead_opt3(mcolor(navy) lcolor(navy) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt3(lcolor(navy%60) lwidth(thin) msize(vsmall)) ///
    lead_ci_opt3(lcolor(navy%60) lwidth(thin) msize(vsmall)) ///
    /// 4) TWFE OLS in emerald
    lag_opt4(mcolor(emerald) lcolor(emerald) msymbol(circle) msize(vsmall)) ///
    lead_opt4(mcolor(emerald) lcolor(emerald) msymbol(circle) msize(vsmall)) ///
    lag_ci_opt4(lcolor(emerald%60) lwidth(thin) msize(vsmall)) ///
    lead_ci_opt4(lcolor(emerald%60) lwidth(thin) msize(vsmall))
* width(2155): match the resolution of the GUI-exported v2 PNGs.
graph export "`outdir'/Human_AMU.png", width(2155) replace
graph export "`outdir'/Human_AMU.pdf", replace

log close
