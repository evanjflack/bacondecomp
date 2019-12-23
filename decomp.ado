*! 1.0.2 15sep2019 Andrew Goodman-Bacon, Thomas Goldring, Austin Nichols
* changes from 1.0.1 31jul2019:
* fix bug pointed out by Nic Duquette, i.e. anyalways not defined to be zero with no always-treated unit in the panel
* and fixed branching with ddetail as reported by David McLaughlin, i.e. main code executes after baconddtiming and exits with error at foreach var of varlist `x' {
* 1.0.1 31jul2019 Andrew Goodman-Bacon, Thomas Goldring, Austin Nichols
* change from 1.0.0 22jul2019: only show graph options when debug option is shown
*! Bacon decomposition of a Diff-in-Diff model per
*! Goodman-Bacon, Andrew. (2019). "Difference-in-Differences with Variation in Treatment Timing." [http://www.nber.org/papers/w25018].
*! This code is licensed under the CC0 1.0 Universal license.
*! The full legal text as well as a human-readable summary can be accessed at http://creativecommons.org/publicdomain/zero/1.0/
prog bacondecomp, eclass
version 10.2
if replay() {
 syntax [anything] [, EForm(string) Level(real 95) * ]
 eret di, eform(`eform') level(`level') nopv noci `options'
 matlist e(sumdd), tw(20) tit(Bacon Decomposition) format(%12.0g) border(all)
}
else {
 syntax [varlist] [if] [in] [aw fw pw iw/] [, Msymbols(string) MColors(string) MSIZes(string) DDLine(string) noLine noNOte cvc within GROpt(string asis) noLEGend noGRaph stub(string) forcebalance(int 2) compare(string) txtformat(string) txtformat2(string) comparenox DDetail debug * ]
 marksample touse
******************************
*parse the varlist into y (outcome), tr (treatment dummy), x (controls)
******************************
 gettoken y rest: varlist
 gettoken tr x: rest
******************************
* Check data xtset and strongly balanced; issue error messages if not xtset or not balanced
******************************
 qui xtset
 loc i=r(panelvar)
 loc t=r(timevar)
 loc b=r(balanced)
 if "`i'"=="" | "`t'"=="" | "`i'"=="." | "`t'"=="." {
  di as err "Both panel and time variables must be set"
  err 459
  }
 qui xtset
 if "`b'"!="strongly balanced" {
  di as err "Panel must be strongly balanced, i.e. each panel must have the same set of time points"
  }
******************************
*if no controls or weights, and ddetail option specified, do local version of ddtiming
******************************
 if ("`x'"=="" & "`exp'"=="") & "`ddetail'"!="" {
  baconddtiming `y' `tr', i(`i') t(`t') ms(`msymbols') mc(`mcolors') msiz(`msizes') `line' ddli(`ddline') stub(`stub') `options'
  }
 if !("`x'"=="" & "`exp'"=="") & "`ddetail'"!="" {
  di as err "ddetail option currently only works without weights and adjustment for control variables"
  err 198
  }
******************************
* add this functionality in next release, also no-control compare option
******************************
 if ("`x'"=="") & "`ddetail'"=="" {
  di as err "with no control variables, must specify the ddetail option currently (only works without weights)"
  err 198
  }

if !("`x'"=="") {
******************************
*Create variables for graph STUB option in syntax; [stub] is prepended to any variables to be saved
******************************
 * (check stub variable does not exist, or create tempvars if stub not specified)
 if "`stub'"=="" {
  tempvar T C S B R2 g cgroup
  }
 else {
  loc T="`stub'T"
  loc C="`stub'C"
  loc S="`stub'S"
  loc B="`stub'B"
  loc R2="`stub'R2"
  loc g="`stub'gp"
  loc cgroup="`stub'cgroup"
  }
 foreach var in T C {
  qui g str1 ``var''=""
  }
 foreach var in S B R2 cgroup {
  qui g byte ``var''=.
  }
 la var `T' "Bacon decomp 2x2 treated group"
 la var `C' "Bacon decomp 2x2 control group"
 la var `S' "Bacon decomp 2x2 weight"
 la var `B' "Bacon decomp 2x2 coefficient"
 la var `R2' "Bacon decomp 2x2 R-squared"
 cap la var `g' "Bacon decomp timing group"
 la var `cgroup' "Bacon decomp aggregation group"

******************************
*Check strongly balanced in nonmissing obs
******************************
 markout `touse' `i' `t' `y' `tr' `x'

 qui xtset
 loc b=r(balanced)
 if "`b'"!="strongly balanced" {
  di as err "Panel must be strongly balanced, i.e. each panel must have the same set of time points"
  di as err "with nonmissing data on panel id, time, and covariates"
  di as txt "(panels are said to be strongly balanced if each panel contains the same time points,"
  di as txt "weakly balanced if each panel contains the same number of observations but not the same"
  di as txt "time points, and unbalanced otherwise); see also " as smcl "{mansection XT xtset}"
  err 459
  }
******************************
* Check treatment weakly increasing in time
* will change after possibly working out
* details of treatment turning on and off in a future paper
******************************
 tempvar negt first last jump
 qui {
  gen `negt'=-`t'
  bys `touse' `i' (`negt'): g `first'=`tr'[_N]
  bys `touse' `i' (`t'): g `last'=`tr'[_N]
  bys `touse' `i' (`t'): g `jump'=`tr'-`tr'[_n-1] if `touse'
  }
 su `last', mean
 if (r(max)!=1) | !inlist(r(min),0,1) {
  di as err "Treatment variable `tr' does not have a maximum of one in last period"
  di as err "or has a minimum not either one or zero in last period"
  err 459
  }
 su `first', mean
 if (r(min)!=0) | !inlist(r(max),0,1) {
  di as err "Treatment variable `tr' does not have a minimum of zero in first period"
  di as err "or has a maximum not either one or zero in first period"
  err 459
  }
 su `jump', mean
 if (r(min)!=0) | (r(max)!=1) {
  di as err "Treatment variable `tr' does not weakly increase (0->1) over time periods"
  err 459
  }
******************************
 * Create variables for onset and timing groups (for dyads), calculations
 * this is like reverse engineering the t* variable from D
******************************
 tempvar t1 onset never always samp
 qui {
  bys `touse' `i' (`t'): g `always'=(`tr'[1]==`tr'[_N]) & `tr'[1]==1 if `touse'
  bys `touse' `i' (`t'): g `never'=(`tr'[1]==`tr'[_N]) & `tr'[1]==0 if `touse'
  bys `touse' `i' (`t'): g `t1'=`t' if `jump'==1 & `touse'
  bys `touse' `i' (`t1'): g `onset'=`t1'[1] if `touse'
******************************
* Just make a consecutive variable for groups
* to make it easier to loop through and map to observation numbers
******************************
  qui egen long `g'=group(`onset'), lab
  su `g', mean
  loc ntimegps=r(max)
  loc alwaysnever=0
******************************
* put always treated units at the end of the list of groups
******************************
  loc anyalways=0
  count if `always'==1
  if r(N)>0 {
    qui replace `g'=`ntimegps'+1 if `always'==1
   la def `g' `=`ntimegps'+1' "Always", modify
   su `g', mean
   loc ntimegps=r(max)
   loc alwaysnever=`alwaysnever'+1
   loc anyalways=1
   }
******************************
* put never treated units at the end of the list of groups (after always)
******************************
  count if `never'==1
  if r(N)>0 {
   qui replace `g'=`ntimegps'+1 if `never'==1
   la def `g' `=`ntimegps'+1' "Never", modify
   su `g', mean
   loc ntimegps=r(max)
   loc alwaysnever=`alwaysnever'+1
   loc anynever=1
   }
  qui compress `g'
  }
 * Mark sample to exclude cases with missing treatment onset
 markout `touse' `g'

 * Check treatment in {0,1} and restrict to balanced panel
 cap assert inlist(`tr',0,1) if `touse'
 if _rc!=0 {
  di as err "Treatment variable `tr' is not either 0 or 1 in all time periods in sample"
  err 459
  }
 qui levelsof `t' if `touse', loc(ts)
 qui foreach i of loc ts {
  tempvar time`i'
  g byte `time`i''=(`t'==`i') if `touse'
  la var `time`i'' "`t'==`i'"
  loc timedummies `timedummies' `time`i''
  }
 if "`exp'"!="" {
  loc wtexp `"[`weight'=`exp']"'
  }
 else loc exp "1"

******************************
* Estimate DiD model
******************************
 qui xtreg `y' `tr' `x' `timedummies' `wtexp'  if `touse', fe `options'
 tempname DDest origb origv
 scalar `DDest'=_b[`tr']
 mat `origb'=e(b)
 mat `origv'=e(V)
 tempvar p d sumwt sumg samp Dtilde dp ptilde pgjtilde
 tempname VD Vb Vdp Rsq finals share BD Bp Beta finals

 * Frisch-Waugh-Lovell Regression
 qui xtreg `tr' `x' `timedummies' `wtexp' if `touse', fe `options'
 qui predict double `p' if e(sample)
 qui predict double `d' if e(sample), ue

******************************
* collapse  x's and p to the group/year level
******************************
 loc xglist
 foreach v in `x' `p' {
  qui {
   tempvar g`v'
   bys `touse' `g' `t':g double `sumg'=sum(`v'*`exp')
   by `touse'  `g' `t':g double `sumwt'=sum(`exp'/(`v'<.))
	by `touse'  `g' `t':g double `g`v''=`sumg'[_N]/`sumwt'[_N] if `touse'
   cap drop `sumg' `sumwt'
   if "`v'"!="`p'" loc xglist `xglist' `g`v''
   }
  }

******************************
* comparisons across all timing groups
* first tell user what is happening in case the regressions take a long time
******************************
loc index 1
di as txt "Computing decomposition across `ntimegps' timing groups"
if `alwaysnever'>0 {
 di as txt "including " cond("`anyalways'"=="1","an always-treated group","") cond(`alwaysnever'==2," and ","") cond("`anynever'"=="1","a never-treated group","")
 }
loc ncompare=`ntimegps'
if `alwaysnever'>1 loc ncompare=`ntimegps'-1

* l is outer loop
* k loops up to l-1
forv it=2/`ntimegps' {
 forv jt=1/`=`it'-1' {
  loc itstring="`:lab `g' `it''"
  loc jtstring="`:lab `g' `jt''"
  qui g byte `samp' =((`g'==`it')|(`g'==`jt'))&`touse'
  su `samp' if `touse' [aw=`exp'], mean
  scalar `share' = r(mean)

******************************
* this is VD so it shouldn't have the x's partialled out
******************************
if !(`alwaysnever'>1 & `it'==`ntimegps' & `jt'==`ntimegps'-1) {
 *get dyad variance
 qui xtreg `tr' `timedummies' `wtexp' if `samp', fe `options'
 qui predict double `Dtilde' if e(sample), e
 qui sum `Dtilde' [aw=`exp']
 scalar `VD' = ((r(N)-1)/r(N))*r(Var)
 }
else {
 scalar `VD' = 0
 }

******************************
* this code partials FE out of the GROUP-level x's, not indiv-level hence `g`var''
******************************
if !(`alwaysnever'>1 & `it'==`ntimegps' & `jt'==`ntimegps'-1) {
 local XXlist
 foreach var of varlist `x' {
  qui xtreg `g`var'' `timedummies' `wtexp' if `samp', fe `options'
  tempvar XX`var'
  qui predict double `XX`var'' if e(sample), e
  local XXlist `XXlist' `XX`var''
  }
 cap drop `pgjtilde'
 qui reg `Dtilde' `XXlist' [aw=`exp']
 scalar `Rsq' = e(r2)
 qui predict double `pgjtilde' if e(sample)
 drop `XXlist'
 *partialled out p in the dyad
 cap drop `ptilde'
 qui xtreg `g`p'' `timedummies' `wtexp' if `samp', fe `options'
 qui predict double `ptilde' if e(sample), e
 qui gen double `dp' = `pgjtilde' - `ptilde' if `samp'
 *get variance of "pg tilde" in the dyad
 qui sum `dp' [aw=`exp']
 scalar `Vdp' = ((r(N)-1)/r(N))*r(Var)
 *weight: this is the variance of "dtilde" in the dyad.
 scalar `finals' = (`share')^2*((1-`Rsq')*`VD' + `Vdp')
 }
else {
 scalar `VD' = 0
 scalar `Vdp' = 0
 scalar `Rsq' = 0
 cap drop `ptilde'
 qui xtreg `g`p'' `timedummies' `wtexp' if `samp', fe `options'
 qui predict double `ptilde' if e(sample), e
 qui sum `ptilde' [aw=`exp']
 scalar `Vb' = ((r(N)-1)/r(N))*r(Var)
 *weight: this is the variance of "dtilde" in the dyad.
 scalar `finals' = (`share')^2*(`Vb')
 }

if !(`alwaysnever'>1 & `it'==`ntimegps' & `jt'==`ntimegps'-1) {
 ****get the proper X-adjusted dyad coef
 tempname BD Bb
 qui xtreg `y' `tr' `xglist' `timedummies' `wtexp' if `samp', fe `options'
 scalar `BD' = _b[`tr']
 qui xtreg `y' `dp' `timedummies' `wtexp' if `samp', fe `options'
 scalar `Bb' = _b[`dp']
 *The dyad "Beta" combines the proper controlled one and a term for how wrong the "FWL" coef is
 scalar `Beta' = ((1-`Rsq')*`VD'*`BD' + `Vdp'*`Bb')/((1-`Rsq')*`VD' + `Vdp')
 }
else {
 ****get the proper X-adjusted dyad coef
 * tempname BD Bb
 qui xtreg `y' `g`p'' `timedummies' `wtexp' if `samp', fe `options'
 scalar `Beta' = -_b[`g`p'']
 }

 qui {
  replace `T' = "`itstring'" in `index'
  replace `C' = "`jtstring'" in `index'
  replace `S' = `finals' in `index'
  replace `B' = `Beta' in `index'
  replace `R2' = `Rsq' in `index'
  replace `cgroup'=1*(!inlist("`itstring'","Always","Never")&!inlist("`jtstring'","Always","Never"))+2*inlist("Always","`itstring'","`jtstring'")+3*inlist("Never","`itstring'","`jtstring'")  in `index'
  loc index=`index'+1
  }
 foreach name in samp Dtilde dp {
   cap drop ``name''
   }
  }
 }
lab def `cgroup' 1 "Timing groups" 2 "Always treated vs timing" 3 "Never treated vs timing" 4 "Within" 5 "Always vs never treated"
la val `cgroup' `cgroup'

******************************
* within part
******************************
tempname Bw Vw
tempvar pw
if `c(matsize)'<11000 cap set matsize 11000
qui xtreg `y' i.`g'##(`timedummies') `p' `wtexp' if `touse', fe `options'
scalar `Bw' = -_b[`p']
qui xtreg `p' i.`g'##(`timedummies') `wtexp' if `touse', fe `options'
qui predict double `pw' if e(sample), e
qui sum `pw' [aw=`exp']
scalar `Vw' = ((r(N)-1)/r(N))*r(Var)
qui replace `T' = "Within" in `index'
qui replace `C' = "" in `index'
qui replace `S' = `Vw' in `index'
qui replace `B' = `Bw' in `index'
qui replace `cgroup'=4 in `index'

 tempvar omega DD wgt sg
 tempname tb totals O
******************************
* rescale weights
******************************
 su `S' if (`T'~="Within"), mean
 scalar `tb'=r(sum)
 qui gen double `sg' = `S'/scalar(`tb') if `T'~="Within"
 su `S', mean
 scalar `totals' = r(sum)
 qui g double `omega'= 1-scalar(`tb')/scalar(`totals')
 su `omega', mean
 scalar `O' = r(sum)
 qui gen double `DD' = `sg'*(1-`O')*`B'*(`T'~="Within")+`O'*`B'*(`T'=="Within") if !mi(`S')
 qui gen double `wgt' = `S'/`totals' if !mi(`S')

******************************
 * Post estimates
******************************
 tempname postb postv diagv output output1 output2 sumcalc summary summary1 summary2 summary3
 forv ix=1/`index' {
  mat `sumcalc'=nullmat(`sumcalc')\(`=`cgroup'[`ix']',`=`B'[`ix']',`=`wgt'[`ix']')
  if (`=`cgroup'[`ix']')==4 loc winest=`=`B'[`ix']'
  if (`=`cgroup'[`ix']')==4 loc winwgt=`=`wgt'[`ix']'
  if (`=`cgroup'[`ix']')==5 loc cvcest=`=`B'[`ix']'
  if (`=`cgroup'[`ix']')==4 loc cvcwgt=`=`wgt'[`ix']'
  mat `postb'=nullmat(`postb'),`=`B'[`ix']'
  mat `postv'=nullmat(`postv'),`=`wgt'[`ix']'
  mat `output1'=nullmat(`output1'),`=`B'[`ix']'
  mat `output2'=nullmat(`output2'),`=`wgt'[`ix']'
  loc rname `rname' "`=`T'[`ix']'`=cond("`=`C'[`ix']'"=="","","_")'`=`C'[`ix']'"
  }
 if "`debug'"!="" di "sumcalc matrix"
 if "`debug'"!="" mat li `sumcalc'
 mata:st_matrix("`summary2'",panelsum(st_matrix("`sumcalc'")[.,3],panelsetup(st_matrix("`sumcalc'"),1)))
 mata:st_matrix("`summary1'",panelsum(st_matrix("`sumcalc'")[.,2],st_matrix("`sumcalc'")[.,3],panelsetup(st_matrix("`sumcalc'"),1)):/panelsum(st_matrix("`sumcalc'")[.,3],panelsetup(st_matrix("`sumcalc'"),1)))
 mata:st_matrix("`summary3'",panelsum(st_matrix("`sumcalc'")[.,1],panelsetup(st_matrix("`sumcalc'"),1)):/panelsum(st_matrix("`sumcalc'")[.,1]:<.,panelsetup(st_matrix("`sumcalc'"),1)))
 mat `summary'=`summary1',`summary2'
 if "`debug'"!=""  di "summary matrix"
 if "`debug'"!=""  mat li `summary'
 forv row=1/`=rowsof(`summary3')' {
  * T vs Never; T vs Always, Timing, Within, C vs C
  if `summary3'[`row',1]==1 loc sumrows `sumrows' "Timing_groups"
  if `summary3'[`row',1]==2 loc sumrows `sumrows' "Always_v_timing"
  if `summary3'[`row',1]==3 loc sumrows `sumrows' "Never_v_timing"
  if `summary3'[`row',1]==4 loc sumrows `sumrows' "Within"
  if `summary3'[`row',1]==5 loc sumrows `sumrows' "Always_v_never"
  }
 mat rownames `summary'=`sumrows'
 mat colnames `summary'=Beta TotalWeight
 mat `output'=`output1'',`output2''
 mat `diagv'=diag(`postv')
 mat rownames `postb'=y1
 mat colnames `postb'=`rname'
 mat rownames `output'=`rname'
 mat colnames `output'=DDest Weight
 mat colnames `postv'=`rname'
 mat rownames `diagv'=`rname'
 mat colnames `diagv'=`rname'
 qui count if `touse'
 loc N = r(N)
 tempname subsetb subsetv
 matrix `subsetb' = `origb'[1,1]
 matrix `subsetv' = `origv'[1,1]
 mat rownames `subsetb'=y1
 mat colnames `subsetb'=`tr'
 mat rownames `subsetv'=`tr'
 mat colnames `subsetv'=`tr'
 if "`debug'"!="" di "subsetb matrix"
 if "`debug'"!=""  mat li `subsetb'
 if "`debug'"!="" di "subsetv matrix"
 if "`debug'"!=""  mat li `subsetv'
 eret post `subsetb' `subsetv', esample(`touse')
 ereturn scalar N = `N'
 ereturn matrix dd = `postb'
 ereturn matrix wt = `postv'
 ereturn matrix sumdd = `summary', copy
 ereturn local depvar "`y'"
 ereturn local version "1.0.0"
 ereturn local cmd "bacondecomp"
 ereturn local properties "b V"
 eret di, eform(`eform') level(`level')
 matlist `summary', tw(20) tit(Bacon Decomposition) format(%12.0g) border(all)
 if "`debug'"!="" di "output matrix"
 if "`debug'"!=""  mat li `output'

******************************
* Graphs, if not suppressed with nograph option
******************************
if "`graph'"=="" {
 loc opt ytitle("2x2 DD Estimate") xtitle("Weight") ylabel(, nogrid) graphregion(fcolor(white) color(white) icolor(white) )
 *plotregion(margin(medsmall)) margin(small)
 if "`line'"=="" loc opt `opt' yline(`=`DDest'', lwidth(thick) `ddline')
 loc capt `" "Overall DD Estimate = `: di `txtformat' `DDest''" "'
 if `alwaysnever'==2 & "`cvc'"=="" {
  loc capt `"`capt' "Always vs never treated = `: di `txtformat' `cvcest'' (weight = `: di `txtformat' `cvcwgt'')""'
  }
 if "`within'"=="" {
  loc capt `"`capt' "Within component = `: di `txtformat' `winest'' (weight = `: di `txtformat' `winwgt'')""'
  }
 loc opt `"`opt' caption(`capt')"'
 if "`debug'"!=""  di `"`opt'"'
 qui levelsof `cgroup', loc(cgroups)
 loc cc 1
 if "`msymbols'"=="" loc msymbols "oh t x d dh dh dh dh"
 if "`mcolors'"=="" loc mcolors "black gray black black black black black"
 if "`msizes'"=="" loc msizes "large large large large large large large large large"
 foreach c of loc cgroups {
  if "`cvc'"=="" & "`cc'"=="5" continue
  if "`within'"=="" & "`cc'"=="4" continue
  loc sym: word `cc' of `msymbols'
  loc col: word `cc' of `mcolors'
  loc siz: word `cc' of `msizes'
  loc sc `sc' sc `B' `wgt' if `cgroup'==`c', msym(`sym') msize(`siz') mcolor(`col')||
  loc cc=`cc'+1
  }
 if "`legend'"=="" & `alwaysnever'==2 {
  loc opt `opt' leg(lab(1 "Timing groups") lab(2 "Always treated vs timing") lab(3 "Never treated vs timing") lab(4 "Always vs never treated") lab(5 "Within"))
  }
 if "`legend'"=="" & `alwaysnever'==1 & `anyalways'==1 {
  loc opt `opt' leg(lab(1 "Timing groups") lab(2 "Always treated vs timing") lab(3 "Within"))
  }
 if "`legend'"=="" & `alwaysnever'==1 & `anynever'==1 {
  loc opt `opt' leg(lab(1 "Timing groups") lab(2 "Never treated vs timing") lab(3 "Within"))
  }
 if "`legend'"=="" & `alwaysnever'==0 {
  loc opt `opt' leg(off)
  }
 `sc', `opt' `gropt'
 }
}
}
end

* adapted from ddtiming; took out save option and added stub option
* version 0.1  11nov2018  Thomas Goldring, thomasgoldring@gmail.com
/* CC0 license information:
To the extent possible under law, the author has dedicated all copyright and related and neighboring rights
to this software to the public domain worldwide. This software is distributed without any warranty.

This code is licensed under the CC0 1.0 Universal license.  The full legal text as well as a
human-readable summary can be accessed at http://creativecommons.org/publicdomain/zero/1.0/
*/

program define baconddtiming, eclass sortpreserve
  version 11.2

  syntax varlist(min=2 max=2 numeric) [if], [i(varname numeric) ///
    t(varname numeric) Msymbols(string) MColors(string) MSIZes(string) ///
    DDLine(string) noLine stub(string) replace *]

  * Perform checks
  if "`i'" == "" {
    di as err `"A panel variable must be specified using option "i({it:panelvar})""'
    exit 198
  }

  if "`t'" == "" {
    di as err `"A time variable must be specified using option "t({it:timevar})""'
    exit 198
  }

  if "`ddline'" != "" & "`line'" != "" {
    di as err `"Cannot specify options "ddline" and "noline" together"'
    exit 198
  }

  if "`replace'" == "" {
    if `"`savegraph'"' != "" {
      if regexm(`"`savegraph'"', "\.[a-zA-Z0-9]+$") confirm new file `"`savegraph'"'
      else confirm new file `"`savegraph'.gph"'
    }
    if `"`savedata'"' != "" {
      confirm new file `"`savedata'.csv"'
      confirm new file `"`savedata'.do"'
    }
  }

  * Mark sample
  marksample touse
  markout `touse' `i' `t'

  * Parse varlist
  tokenize `varlist'
  local y `1'
  local tr `2'

  * Check treatment variable is binary
  if !inlist(`tr',0,1) {
    di as error "Treatment variable must be binary 0/1"
    exit 198
  }

  * Create temporary names and variables
  tempname untr_samp_share contr_samp_share tr_samp_share total_untr_samp_share ///
    tr_grps tr_share sigma_term1 sigma_term2 sigma ///
    tmp tmp_s tmp_mu tmp_s_mu tmp_s_1mu ///
    wt_untr wt_untr_tot wt_contr wt_contr_tot wt_tr_e wt_tr_e_tot ///
    wt_tr_l wt_tr_l_tot lrange urange savedatafile ///
    dd_untr_est dd_contr_est dd_e dd_l dd_e_tr_est dd_l_tr_est ///
    dd_est dd_wt dd_mat ///
    untr_est_avg contr_est_avg e_tr_est_avg l_tr_est_avg
  tempvar tr_time tr_never tr_before gr_dd_est gr_dd_wt gr_dd_type

  quietly {

  * Calculate treatment time
  noisily di as txt "Calculating treatment times..."
  bys `i' (`t'): gen `tr_time' = `t' * (sum(`tr') == 1 & sum(`tr'[_n-1]) == 0) if `touse'
  bys `i' (`tr_time'): replace `tr_time' = `tr_time'[_N] if `touse'

  * Indicators for never treated and already treated groups
  bys `i' (`t'): gen `tr_never' = (`tr'[_N] == 0 & `touse')
  bys `i' (`t'): gen `tr_before' = (`tr'[1] == 1 & `touse')
  replace `tr_time' = 0 if `tr_before' == 1 & `touse'

  * Indicator for presence of untreated group
  count if `tr_never' == 1
  if r(N) > 0 local untr 1
  else local untr 0

  * Indicator for presence of already treated group
  count if `tr_before' == 1
  if r(N) > 0 local contr 1
  else local contr 0

  * CREATE WEIGHTS

  * Group fraction of sample
  tab `tr_time' if `tr_time' == 0 & `tr_before' == 0 & `touse', matcell(`untr_samp_share')
  tab `tr_time' if `tr_time' == 0 & `tr_before' == 1 & `touse', matcell(`contr_samp_share')
  tab `tr_time' if `tr_time' != 0 & `touse', matcell(`tr_samp_share')

  tab `tr_time' if `touse'
  cap mat `untr_samp_share' = `untr_samp_share' / r(N)
  cap mat `contr_samp_share' = `contr_samp_share' / r(N)
  cap mat `tr_samp_share' = `tr_samp_share' / r(N)

  if `untr' == 1 & `contr' == 1 {
    mat `total_untr_samp_share' = `untr_samp_share' + `contr_samp_share'
  }
  else if `untr' == 1 & `contr' == 0 {
    mat `total_untr_samp_share' = `untr_samp_share'
  }
  else if `untr' == 0 & `contr' == 1 {
    mat `total_untr_samp_share' = `contr_samp_share'
  }

  * Fraction of time spent treated
  tab `tr_time' if `tr_time' != 0 & `touse', matrow(`tr_grps') // Treatment groups by treatment time

  forval k = 1/`= rowsof(`tr_grps')' {
    sum `tr' if `tr_time' == `tr_grps'[`k',1] & `touse', meanonly
    mat `tr_share' = nullmat(`tr_share')\r(mean)
  }

  * Variance of treatment
  if `untr' == 1 | `contr' == 1 {
    forval k = 1/`= rowsof(`tr_grps')' {
      scalar `tmp' = `tr_samp_share'[`k',1] * `total_untr_samp_share'[1,1] * `tr_share'[`k',1] * (1 - `tr_share'[`k',1])
      mat `sigma_term1' = nullmat(`sigma_term1')\\`tmp'
    }
  }
  forval k = 1/`= rowsof(`tr_grps')' {
    local l = `k' + 1
    while `l' <= `= rowsof(`tr_grps')' {
	  mat `tmp' = `tr_samp_share'[`k',1] * `tr_samp_share'[`l',1] * (`tr_share'[`k',1] - `tr_share'[`l',1]) * (1 - (`tr_share'[`k',1] - `tr_share'[`l',1]))
      mat `sigma_term2' = nullmat(`sigma_term2')\\`tmp'
      local ++l
    }
  }

  if `untr' == 1 | `contr' == 1 {
    mata : st_matrix("`sigma'", colsum(st_matrix("`sigma_term1'")) + colsum(st_matrix("`sigma_term2'")))
  }
  else {
    mata : st_matrix("`sigma'", colsum(st_matrix("`sigma_term2'")))
  }

  noisily di as txt "Calculating weights..."

  * Weights where comparison group is untreated (never treated)
  if `untr' == 1 {
    forvalues k = 1/`= rowsof(`tr_grps')' {
      mat `tmp' = `k',0,((`tr_samp_share'[`k',1] * `untr_samp_share'[1,1] * `tr_share'[`k',1] * (1 - `tr_share'[`k',1])) / `sigma'[1,1]),3
      mat `wt_untr' = nullmat(`wt_untr')\\`tmp'
    }
	mata : st_matrix("`wt_untr_tot'", colsum(st_matrix("`wt_untr'")[1...,3]))
	mata : st_matrix("`tmp'", st_matrix("`wt_untr'")[1...,3] / st_matrix("`wt_untr_tot'"))
	mat `wt_untr' = `wt_untr',`tmp'  // Columns: treatment group, control group, weight, comparison type, weight scaled to 1
  }

  * Weights where comparison group is untreated (already treated)
  if `contr' == 1 {
    forvalues k = 1/`= rowsof(`tr_grps')' {
      mat `tmp' = `k',0,((`tr_samp_share'[`k',1] * `contr_samp_share'[1,1] * `tr_share'[`k',1] * (1 - `tr_share'[`k',1])) / `sigma'[1,1]),4
      mat `wt_contr' = nullmat(`wt_contr')\\`tmp'
    }
	mata : st_matrix("`wt_contr_tot'", colsum(st_matrix("`wt_contr'")[1...,3]))
	mata : st_matrix("`tmp'", st_matrix("`wt_contr'")[1...,3] / st_matrix("`wt_contr_tot'"))
	mat `wt_contr' = `wt_contr',`tmp'
  }

  * Weights where comparison group is ever treated
  forvalues k = 1/`= rowsof(`tr_grps')' {
    local l = `k' + 1
    while `l' <= `= rowsof(`tr_grps')' {
      mat `tmp_s' = (`tr_samp_share'[`k',1] * `tr_samp_share'[`l',1] * (`tr_share'[`k',1] - `tr_share'[`l',1]) * (1 - (`tr_share'[`k',1] - `tr_share'[`l',1]))) / `sigma'[1,1]
      mat `tmp_mu' = (1 - `tr_share'[`k',1]) / (1 - (`tr_share'[`k',1] - `tr_share'[`l',1]))

      mat `tmp_s_mu' = `k',`l',(`tmp_s' * `tmp_mu'),1
      mat `tmp_s_1mu' = `l',`k',(`tmp_s' * (1 - `tmp_mu')),2

      mat `wt_tr_e' = nullmat(`wt_tr_e')\\`tmp_s_mu'   // e = Treatment's treatment earlier than control's treatment
      mat `wt_tr_l' = nullmat(`wt_tr_l')\\`tmp_s_1mu'  // l = Treatment's treatment later than control's treatment

      local ++l
    }
  }
  mata : st_matrix("`wt_tr_e'", sort(st_matrix("`wt_tr_e'"), (1,2)))
  mata : st_matrix("`wt_tr_l'", sort(st_matrix("`wt_tr_l'"), (1,2)))

  mata : st_matrix("`wt_tr_e_tot'", colsum(st_matrix("`wt_tr_e'")[1...,3]))
  mata : st_matrix("`tmp'", st_matrix("`wt_tr_e'")[1...,3] / st_matrix("`wt_tr_e_tot'"))
  mat `wt_tr_e' = `wt_tr_e',`tmp'  // Columns: treatment group, control group, weight, comparison type, weight scaled to 1

  mata : st_matrix("`wt_tr_l_tot'", colsum(st_matrix("`wt_tr_l'")[1...,3]))
  mata : st_matrix("`tmp'", st_matrix("`wt_tr_l'")[1...,3] / st_matrix("`wt_tr_l_tot'"))
  mat `wt_tr_l' = `wt_tr_l',`tmp'

  * ESTIMATE 2x2 DIFF-IN-DIFF REGRESSIONS

  noisily di as text "Estimating 2x2 diff-in-diff regressions..."

  * Never treated diff-in-diff estimates
  if `untr' == 1 {
    forvalues k = 1/`= rowsof(`tr_grps')' {
      sum `t'
      mat `tmp' = `k',0,`tr_grps'[`k',1],0,r(min),r(max),3
      areg `y' `tr' i.`i' if (`tr_never' == 1 | `tr_time' == `tr_grps'[`k',1]) & `touse', a(`t')
      mat `tmp' = `tmp',_b[`tr']
      mat `dd_untr_est' = nullmat(`dd_untr_est')\\`tmp'
     // Columns: (1) treatment group (2) control group (3) treatment's treatment time
     //  (cont.) (4) control's treatment time (5) start time (6) end time
     //  (cont.) (7) comparison category (8) diff-in-diff estimate
    }
  }

  * Already treated diff-in-diff estimates
  if `contr' == 1 {
    forvalues k = 1/`= rowsof(`tr_grps')' {
      sum `t'
      mat `tmp' = `k',0,`tr_grps'[`k',1],0,r(min),r(max),4
      areg `y' `tr' i.`i' if (`tr_before' == 1 | `tr_time' == `tr_grps'[`k',1]) & `touse', a(`t')
      mat `tmp' = `tmp',_b[`tr']
      mat `dd_contr_est' = nullmat(`dd_contr_est')\\`tmp'
    }
  }

  * Calculate lower and upper ranges for treatment time (only one range is used in each regression)
  sum `t'
  forvalues k = 1/`= rowsof(`tr_grps')' {
    mat `lrange' = nullmat(`lrange')\[r(min),`tr_grps'[`k',1] - 1]
    mat `urange' = nullmat(`urange')\[`tr_grps'[`k',1],r(max)]
    // Columns: (1) time lower bound (2) time upper bound
  }
  forvalues k = 1/`= rowsof(`tr_grps')' {           // Treatment index
    forvalues l = 1/`= rowsof(`tr_grps')' {         // Control index
      if `k' != `l' {
        mat `tmp' = `k',`l',`tr_grps'[`k',1],`tr_grps'[`l',1]

        if `tmp'[1,3] < `tmp'[1,4] {                // If treatment's treatment time < control's treatment time
          mat `tmp' = `tmp',`lrange'[`l',1...],1    // ... use range below control's treatment time
        }
        else mat `tmp' = `tmp',`urange'[`l',1...],2 // ... else use range above control's treatment time

        if `k' < `l' mat `dd_e' = nullmat(`dd_e')\\`tmp'
        else mat `dd_l' = nullmat(`dd_l')\\`tmp'
      }
    }
  }

  local earlylate "e l"
  foreach x of local earlylate {
    forvalues k = 1/`= rowsof(`dd_`x'')' {
      areg `y' `tr' i.`i' if ///
        (`tr_time' == `dd_`x''[`k',3] | `tr_time' == `dd_`x''[`k',4]) & ///
        `t' >= `dd_`x''[`k',5] & `t' <= `dd_`x''[`k',6] & `touse', a(`t')
      mat `dd_`x'_tr_est' = nullmat(`dd_`x'_tr_est')\_b[`tr']
    }
    mat `dd_`x'_tr_est' = `dd_`x'',`dd_`x'_tr_est'
    mata : st_matrix("`dd_`x'_tr_est'", sort(st_matrix("`dd_`x'_tr_est'"), (1,2)))
  }

  * Stitch DD and weights matrices together
  mat `dd_est' = `dd_e_tr_est'\\`dd_l_tr_est'
  mat `dd_wt' = `wt_tr_e'\\`wt_tr_l'

  if `untr' == 1 {
    mat `dd_est' = `dd_est'\\`dd_untr_est'
    mat `dd_wt' = `dd_wt'\\`wt_untr'
  }
  if `contr' == 1 {
    mat `dd_est' = `dd_est'\\`dd_contr_est' // 2x2 DD estimates matrix
    mat `dd_wt' = `dd_wt'\\`wt_contr'       // Weights matrix
  }

  * Put matrices in tempvars for graphing
  gen `gr_dd_est' = .
  gen `gr_dd_wt' = .
  gen `gr_dd_type' = .

  forval k = 1/`= rowsof(`dd_est')' {
    replace `gr_dd_est' = `dd_est'[`k',8] in `k'
    replace `gr_dd_wt' = `dd_wt'[`k',3] in `k'
    replace `gr_dd_type' = `dd_est'[`k',7] in `k'
  }

  } // Ends quietly

  * Calculate DD estimate
  mata : st_matrix("`dd_mat'", st_matrix("`dd_est'")[1...,8] :* st_matrix("`dd_wt'")[1...,3])
  mata : st_matrix("`dd_mat'", colsum(st_matrix("`dd_mat'"))) // Two-way FE estimate

  * Calculate weighted DD estimates by comparison type
  if `untr' == 1 {
    mata : st_matrix("`untr_est_avg'", st_matrix("`dd_untr_est'")[1...,8] :* st_matrix("`wt_untr'")[1...,5])
    mata : st_matrix("`untr_est_avg'", colsum(st_matrix("`untr_est_avg'")))
  }

  if `contr' == 1 {
    mata : st_matrix("`contr_est_avg'", st_matrix("`dd_contr_est'")[1...,8] :* st_matrix("`wt_contr'")[1...,5])
    mata : st_matrix("`contr_est_avg'", colsum(st_matrix("`contr_est_avg'")))
  }

  mata : st_matrix("`e_tr_est_avg'", st_matrix("`dd_e_tr_est'")[1...,8] :* st_matrix("`wt_tr_e'")[1...,5])
  mata : st_matrix("`l_tr_est_avg'", st_matrix("`dd_l_tr_est'")[1...,8] :* st_matrix("`wt_tr_l'")[1...,5])

  mata : st_matrix("`e_tr_est_avg'", colsum(st_matrix("`e_tr_est_avg'")))
  mata : st_matrix("`l_tr_est_avg'", colsum(st_matrix("`l_tr_est_avg'")))

  * PREPARE GRAPH

  * Define scatter command
  forval k = 1/4 {
    local scatter`k' "scatter `gr_dd_est' `gr_dd_wt' if `gr_dd_type' == `k'"
  }

  * Define labels
  local legend1 `"lab(1 "Earlier Group Treatment vs. Later Group Control")"'
  local legend2 `"lab(2 "Later Group Treatment vs. Earlier Group Control")"'

  if `untr' == 1 local legend3 `"lab(3 "Treatment vs. Never Treated")"'

  if `untr' == 1 & `contr' == 1 local legend4 `"lab(4 "Treatment vs. Already Treated")"'
  if `untr' == 0 & `contr' == 1 local legend4 `"lab(3 "Treatment vs. Already Treated")"'

  * Define marker symbols
  local msym1 "X"
  local msym2 "X"
  local msym3 "T"
  local msym4 "Oh"

  if "`msymbols'" != "" {
    forval k = 1/`: word count `msymbols'' {
      local msym`k' : word `k' of `msymbols'
    }
  }

  * Define marker colors
  local mcol1 "gs8"
  local mcol2 "black"
  local mcol3 "gs6"
  local mcol4 "black"

  if "`mcolors'" != "" {
    forval k = 1/`: word count `mcolors'' {
      local mcol`k' : word `k' of `mcolors'
    }
  }

  * Define marker sizes
  local msiz1 "medium"
  local msiz2 "medium"
  local msiz3 "medium"
  local msiz4 "medium"

  if "`msizes'" != "" {
    forval k = 1/`: word count `msizes'' {
      local msiz`k' : word `k' of `msizes'
    }
  }

  forval k = 1/4 {
    local scatter`k' "(`scatter`k'', ms(`msym`k'') mc(`mcol`k'') msiz(`msiz`k''))"
  }
  if `untr' == 0 local scatter3 ""
  if `contr' == 0 local scatter4 ""

  if "`line'" == "" {
    local gr_dd = `dd_mat'[1,1]
    local yline "yline(`gr_dd',`ddline')"
  }
  else local yline ""

  local graphcmd tw `scatter1' `scatter2' `scatter3' `scatter4', ///
    xlabel(,format(%5.2f)) ytitle("2x2 DD Estimate") xtitle("Weight") ///
	`yline' graphregion(color(white)) legend(col(1) `legend_order' ///
	`legend1' `legend2' `legend3' `legend4') `options'
  `graphcmd'

  local graphsave tw (scatter dd_est dd_wt if dd_type == 1) ///
    (scatter dd_est dd_wt if dd_type == 2)

  if `untr' == 1 local graphsave `graphsave' (scatter dd_est dd_wt if dd_type == 3)
  if `contr' == 1 local graphsave `graphsave' (scatter dd_est dd_wt if dd_type == 4)

  local graphsave_opts xlabel(,format(%5.2f)) ytitle("2x2 DD Estimate") xtitle("Weight") ///
    `yline' graphregion(color(white)) legend(col(1) `legend_order' ///
	`legend1' `legend2' `legend3' `legend4') `options'

  * Save graph
  if `"`savegraph'"' != "" {
    if regexm(`"`savegraph'"',"\.[a-zA-Z0-9]+$") local graphextension = regexs(0) /// Check file extension using a regular expression
    if inlist(`"`graphextension'"',".gph","") graph save `"`savegraph'"', `replace'
    else graph export `"`savegraph'"', `replace'
  }

  * Return scalars
  ereturn clear
  ereturn scalar dd = `dd_mat'[1,1]
  ereturn scalar dd_avg_e = `e_tr_est_avg'[1,1]
  ereturn scalar dd_avg_l = `l_tr_est_avg'[1,1]
  if `untr' == 1 {
    ereturn scalar dd_avg_u = `untr_est_avg'[1,1]
  }
  if `contr' == 1 {
    ereturn scalar dd_avg_a = `contr_est_avg'[1,1]
  }
  ereturn scalar wt_sum_e = `wt_tr_e_tot'[1,1]
  ereturn scalar wt_sum_l = `wt_tr_l_tot'[1,1]
  if `untr' == 1 {
    ereturn scalar wt_sum_u = `wt_untr_tot'[1,1]
  }
  if `contr' == 1 {
    ereturn scalar wt_sum_a = `wt_contr_tot'[1,1]
  }

  * Print output
  di as smcl as txt ""
  di as smcl as txt "Diff-in-diff estimate: " as res %-9.3f `dd_mat'[1,1]
  di as smcl as txt ""
  di as smcl as txt "DD Comparison              Weight      Avg DD Est"
  di as smcl as txt "{hline 49}"
  di as smcl as txt "Earlier T vs. Later C       " as res %5.3f `wt_tr_e_tot'[1,1] "       " as res %9.3f `e_tr_est_avg'[1,1]
  di as smcl as txt "Later T vs. Earlier C       " as res %5.3f `wt_tr_l_tot'[1,1] "       " as res %9.3f `l_tr_est_avg'[1,1]
  if `untr' == 1 {
    di as smcl as txt "T vs. Never treated         " as res %5.3f `wt_untr_tot'[1,1] "       " as res %9.3f `untr_est_avg'[1,1]
  }
  if `contr' == 1 {
    di as smcl as txt "T vs. Already treated       " as res %5.3f `wt_contr_tot'[1,1] "       " as res %9.3f `contr_est_avg'[1,1]
  }
  di as smcl as txt "{hline 49}"
  di as smcl as txt "T = Treatment; C = Control"

  * Save data
  if "`stub'" != "" {
    * "dd_est,weight,weight_rescale,time_lower,time_upper,dd_type" columns
    foreach v in dd_est,weight,weight_rescale,time_lower,time_upper,dd_type {
     g `stub'`v'=.
     }
    forvalues k = 1/`= rowsof(`dd_est')' {
     replace `stub'dd_est = (`dd_est'[`k',8]) in `k'
     replace `stub'weight = (`dd_wt'[`k',3]) in `k'
     replace `stub'weight_rescale = (`dd_wt'[`k',5])
     replace `stub'time_lower = (`dd_est'[`k',5])
     replace `stub'time_upper =  (`dd_est'[`k',6])
     replace `stub'dd_type =  (`dd_est'[`k',7])
     }
    label var `stub'dd_est "2x2 DD estimate"
    label var `stub'weight "2x2 DD weight"
    label var `stub'weight_rescale "2x2 DD weight rescaled within comparison type"
    label var `stub'time_lower "2x2 DD time lower bound"
    label var `stub'time_upper "2x2 DD time upper bound"
    label var `stub'dd_type "2x2 DD type"
    label define `stub'dd_type 1 "Earlier T vs. Later C" 2 "Later T vs. Earlier C" 3 "T vs. Never T" 4 "T vs. Already T"
    label values `stub'dd_type `stub'dd_type
    }

  * Save data
  if "`savedata'" != "" {
    file open `savedatafile' using `"`savedata'.csv"', write text `replace'
    file write `savedatafile' "dd_est,weight,weight_rescale,time_lower,time_upper,dd_type" _n // Row headers

    forvalues k = 1/`= rowsof(`dd_est')' {
      file write `savedatafile' (`dd_est'[`k',8]) "," (`dd_wt'[`k',3]) "," (`dd_wt'[`k',5]) "," (`dd_est'[`k',5]) "," (`dd_est'[`k',6]) "," (`dd_est'[`k',7]) _n
    }

    file close `savedatafile'
    di as smcl as txt ""
    di as smcl as txt `"File `savedata'.csv written containing saved data"'

    * Save a do-file with the commands to generate a labeled dataset and re-create the ddtiming graph
    file open `savedatafile' using `"`savedata'.do"', write text `replace'
    file write `savedatafile' `"insheet using `savedata'.csv"' _n _n

    file write `savedatafile' `"label var dd_est "2x2 DD estimate""' _n
    file write `savedatafile' `"label var weight "2x2 DD weight""' _n
	file write `savedatafile' `"label var weight_rescale "2x2 DD weight rescaled within comparison type""' _n
    file write `savedatafile' `"label var time_lower "2x2 DD time lower bound""' _n
    file write `savedatafile' `"label var time_upper "2x2 DD time upper bound""' _n
    file write `savedatafile' `"label var dd_type "2x2 DD types""' _n _n
    file write `savedatafile' `"label define dd_type 1 "Earlier T vs. Later C" 2 "Later T vs. Earlier C" 3 "T vs. Never T" 4 "T vs. Already T""' _n
    file write `savedatafile' `"label values dd_type dd_type"' _n _n

	file write `savedatafile' `"`graphsave'"' `"`graphsave_opts'"' _n

    file close `savedatafile'
    di as smcl as txt `"File `savedata'.do written containing commands to process saved data"'

  }

end

exit
