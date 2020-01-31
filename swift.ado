/*==================================================
project:       SWIFT (Survey of Well-being via Instant and Frequent Tracking)
Author:        R.Andres Castaneda 
E-email:       acastanedaa@worldbank.org
url:           
Dependencies:  The World Bank
----------------------------------------------------
Creation Date:    30 Jan 2020 - 14:29:44
Modification Date:   
Do-file version:    01
References:          
Output:             
==================================================*/

/*==================================================
0: Program set up
==================================================*/
program define swift, rclass
syntax varname [aw pw iw], ///
POVLine(numlist)           ///
Xvars(varlist)             ///
[                          ///
Nsplit(integer 10)         ///
CLUster(varname)           ///
stpe(real 0.005)           /// Starting pe in stepwise
LIMITpe(real 0.1)          /// Max pe allowed
DELTApe(real 0.005)        /// delta in pe
TOLERANCEpe(real 0.00001)  /// Tolenrance of pr to pe (local pr = `pe' + `tolerancepe') 
ADDm(integer 20)           /// number of m in mi 
seed(integer 12345678)     /// seed
noCHART                    ///
GENWelfare(string)         ///
GENPov(string)             ///
GENFold(string)            ///
SAVEGraph(string)          ///
]
version 16

//------------Conditions
marksample touse 

qui {
  
  cap which midiagplots
  if (_rc) net install st0263.pkg
  
  
  /*==================================================
  1: preapre
  ==================================================*/
  
  //------------Weigts
  tempvar wi 
  if missing(`"`weight'"') generate byte `wi' = 1
  else generate `wi' `exp'
  local wgtcall "[`weight' = `wi']"
  
  //------------SVY set
  svyset `wgtcall'
  
  *----------1.1: Split sample
  tempvar fold
  splitsample `varlist', gen(`fold') nsplit(`nsplit') rround cluster(`cluster')
  
  *----------1.2: Real poverty
  tempvar pov
  gen `pov' = (`varlist' < `povline') if `touse'
  
  
  /*==================================================
  2: Loop 
  ==================================================*/
  tempfile basef
  save `basef'
  
  tempname X
  
  noi disp _n in y "Progress: "
  
  foreach i of numlist 1/`nsplit' {
    
    noi disp _n in y "Fold `i' " _c
    // sample conditions
    tempvar bifold 
    gen `bifold' = (`fold' ==`i')
    tempfile baset
    save `baset'
    
    
    // Actual poverty rate in fold i
    svy: mean `pov' `iff_ts'
    local rpov = e(b)[1,1]	
    
    
    local pe = `stpe'
    local pr = `pe' + `tolerancepe'
    
    while `pe' <= `limitpe' {
      noi disp in w "." _c
      
      // stepwise, get variables
      cap {
        swift_estimate `varlist' `wgtcall', povline(`povline') xvars(`xvars') /* 
        */ pe(`pe') pr(`pr') addm(`addm') seed(`seed') bifold(`bifold') 
        
        local mipov = `r(mipov)'
        
        local diffmse = abs(`rpov' - `mipov')
        
        matrix `X' = nullmat(`X') \ `pe', `i', `rpov', `r(r2_a)', `r(mse)', `mipov', `diffmse'
      }
      if (_rc) {
        noi disp in r "error occurred during fold `i' and pe `pe'" _n 
      }
      
      local pe = `pe' + `deltape'
      local pr = `pe' + `tolerancepe'
      
      use `baset', clear
      
    }
    use `basef', clear
    
  }
  disp _n ""
  /*==================================================
  3: Organize results
  ==================================================*/
  
  //------------return matrix X
  mat colnames `X' = pe fold_ poor r2 mse imp_poor absdiff
  return matrix pmi = `X', copy
  
  //------------calculate mean of MSE  and diff in poverty rates by pe
  drop _all 
  svmat `X', n(col)
  gen n = _n
  
  // P where the absolute value of the difference between the actual and 
  //the projected poverty rates is minimized
  sum absdiff, meanonly
  scalar mindiff = r(min)  // minimum diff in poverty rates
  
  sum n if absdiff == mindiff, meanonly
  local n = r(max)
  local mindiff_pe = `: disp pe[`n']'
  return local mindiff_pe = `mindiff_pe'
  
  // P where the meas squared error is minimized
  sum mse, meanonly
  scalar minmse = r(min)  // minimum diff in poverty rates
  
  sum n if mse == minmse, meanonly
  local n = r(max)
  local minmse_pe = `: disp pe[`n']'
  return local minmse_pe = `minmse_pe'
  
  collapse (mean) mmse = mse (mean) mabsdiff = absdiff, by(pe)
  
  tempname Y
  mkmat *, matrix(`Y')
  
  return matrix mmse = `Y'
  
  if ("`chart'" != "nochart") {
    
    if ("`savegraph'" == "") {
      local savegraph "cross_validation.gph"
    } 
    else {
      if !regexm(lower("`savegraph'"), "\.[a-z]+$")  {
        local savegraph  = "`savegraph'" + ".gph"
      }
    }
    
    tempname msec absdiffc 
    twoway scatter mmse pe, name(`msec', replace)
    twoway scatter mabsdiff pe, name(`absdiffc', replace)
    gr combine `msec' `absdiffc', name(cross_validation, replace)
    gr save "`savegraph'", replace
  }
  
  
  //========================================================
  // predict povert in the whole sample
  //========================================================
  
  use `basef', clear
  
  if ("`genfold'" == "") {
    local genfold "bifold"
  }
  
  expand 2, gen(`genfold')
  
  if ("`genwelfare'" == "") {
    local genwelfare "`varlist'_mi"
  }
  if ("`genpov'" == "") {
    local genpov "pov_mi"
  }
  
  local pr = `mindiff_pe' + `tolerancepe'
  
  swift_estimate `varlist' `wgtcall', povline(`povline') xvars(`xvars') /* 
  */ pe(`mindiff_pe') pr(`pr') addm(`addm') seed(`seed') /* 
  */ bifold(`genfold')  welfimp(`genwelfare') povimp(`genpov')
  return add
  
  return local cmd_kpob = "keep if (_mi_m != 0 | `genfold' == 1)"
  return local cmd_mean_aft = "mi estimate: mean `genpov' `wgtcall'"
  return local cmd_mean_bfr = "mi estimate: mean `genpov' `wgtcall' if (_mi_m != 0 | `genfold' == 1)"
  
} // end of qui 


end
exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

Notes:
1.The optimal p-value is the value where the absolute value of the difference between the actual and the projected poverty rates is minimized. The mean squared error is also examined to check whether the over-fitting problem occurs. If the mean squared error is minimized at a level of p that is smaller than the value where the absolute difference between the actual and the projected poverty rates is minimized, then the former value is chosen as the optimal number.
2.
3.


Version Control:





