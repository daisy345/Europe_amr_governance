* AAMR staggered-DID event study, K=-10 through K=14.
* Third palette trial: named colors red / orange / navy / emerald (CS / Gardner / Stacked / TWFE).

capture log close _all
clear all
version 18.0
set more off
local projectdir "`c(pwd)'"
adopath ++ "`projectdir'/code/vendor/csdid"
adopath ++ "`projectdir'/code/vendor/drdid"
local datafile "`projectdir'/input/Animal_AMR_score.dta"
local outdir "`projectdir'/output"
local kmin = -10
local kmax = 14
log using "`outdir'/Animal_AMR_score.log", replace text

foreach cmd in ftools reghdfe did2s event_plot csdid stackedev {
    cap which `cmd'
    if _rc {
        di as error "Required command `cmd' is unavailable."
        exit 199
    }
}

use "`datafile'", clear
rename amr Y
tsset i t
by i (t): egen Ei = min(cond(D == 1, t, .))
gen K = t - Ei

* Event-time dummies within the requested display/estimation window;
* K=-1 is the omitted reference period.
forvalues k = `kmin'/-2 {
    local a = abs(`k')
    gen F`a'event = (K == `k')
}
forvalues k = 0/`kmax' {
    gen L`k'event = (K == `k')
}

* 1. Callaway and Sant'Anna (2020)
gen gvar = cond(missing(Ei), 0, Ei)
csdid Y, ivar(i) time(t) gvar(gvar) agg(event)
matrix cs_b = e(b)
matrix cs_v = e(V)

* 2. Gardner DID-2S (2021)
did2s Y, first_stage(i.i i.t) second_stage(F*event L*event) treatment(D) cluster(i)
matrix did2s_b = e(b)
matrix did2s_v = e(V)

* 3. Stacked event-study DID (Cengiz et al., 2019)
gen treat_year = Ei
gen no_treat = missing(Ei)
foreach var of varlist F*event L*event {
    replace `var' = 0 if no_treat
}
preserve
stackedev Y F*event L*event, cohort(treat_year) time(t) ///
    never_treat(no_treat) unit_fe(i) clust_unit(i)
restore
matrix stacked_b = e(b)
matrix stacked_v = e(V)

* 4. Conventional TWFE event study
reghdfe Y F*event L*event, absorb(i t) vce(cluster i)
estimates store ols

* Combined figure: third palette trial (red/orange/navy/emerald), vsmall
* filled-circle markers, thin 60%-opacity CI lines with short caps,
* series dodged by -0.21(0.08)0.21.
event_plot cs_b#cs_v did2s_b#did2s_v stacked_b#stacked_v ols, ///
    stub_lag(Tp# L#event L#event L#event) ///
    stub_lead(Tm# F#event F#event F#event) ///
    plottype(scatter) ciplottype(rcap) together ///
    perturb(-0.21(0.08)0.21) noautolegend trimlead(10) trimlag(14) ///
    graph_opt( ///
    xtitle("Time since NAP launch (years)") ///
    ytitle("Animal AMR score") xlabel(`kmin'(1)`kmax', nogrid) ///
    xscale(range(`kmin' `kmax') extend) ylabel(, angle(horizontal) nogrid) ///
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
graph export "`outdir'/Animal_AMR_score.png", width(2155) replace
graph export "`outdir'/Animal_AMR_score.pdf", replace

log close
