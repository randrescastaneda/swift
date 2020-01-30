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
]
version 16

//------------Conditions
preserve
marksample touse 

qui {
  
  
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
    local iff_tr "if (`touse' & `fold' !=`i')" // training
    local iff_ts "if (`touse' & `fold' ==`i')" // test
    
    // Actual poverty rate in fold i
    svy: mean `pov' `iff_ts'
    local rpov = e(b)[1,1]	
    
    
    local pe = `stpe'
    local pr = `pe' + `tolerancepe'
    
    while `pe' <= `limitpe' {
      noi disp in w "." _c
      
      // stepwise, get variables
      cap {
        stepwise, pr(`pr') pe(`pe'): reg `varlist' `xvars' `wgtcall' `iff_tr'
        local xvarm: colnames e(b)
        local xvarm: subinstr local xvarm "_cons" "", word
        
        // regression in training
        reg `varlist' `xvarm' `wgtcall' `iff_tr'
        local r2_a = e(r2_a)
        
        // Predict in test
        tempvar yhat sqerr
        predict `yhat' `iff_ts'
        gen `sqerr' = (`varlist' - `yhat')^2 `iff_ts'
        svy: mean `sqerr' `iff_ts'
        local mse = e(b)[1,1]
        
        // set to missing all welfare values for mi
        replace `varlist' =. `iff_ts'
        
        // Multiple imputation
        mi set mlong                       // define output
        mi register imputed `varlist'      // register welfare variabe
        mi impute reg `varlist' `xvarm' `wgtcall', add(`addm') rseed(`seed') force
        
        // new variabe
        tempvar imp_w
        gen `imp_w' = (`varlist' < `povline' ) & _mi_m > 0
        mi estimate, post: svy: mean `imp_w' `iff_ts'
        local mipov = e(b)[1,1]
        
        local diffmse = abs(`rpov' - `mipov')
        
        matrix `X' = nullmat(`X') \ `pe', `i', `rpov', `r2_a', `mse', `mipov', `diffmse'
      }
      if (_rc) {
        noi disp in r "error occurred during folde `i' and pe `pe'" _n 
      }
      
      local pe = `pe' + `deltape'
      local pr = `pe' + `tolerancepe'
      
      use `basef', clear
    }
    
  }
  
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
  local mindiff_pe = r(max)
  return local mindiff_pe = `mindiff_pe'
  
  // P where the meas squared error is minimized
  sum mse, meanonly
  scalar minmse = r(min)  // minimum diff in poverty rates
  
  sum n if mse == minmse, meanonly
  local minmse_pe = r(max)
  return local minmse_pe = `minmse_pe'
  
  collapse (mean) mmse = mse (mean) mabsdiff = absdiff, by(pe)
  
  tempname Y
  mkmat *, matrix(`Y')
  
  return matrix mmse = `Y'
  
  if ("`chart'" != "nochart") {
    tempname msec absdiffc 
    twoway scatter mmse pe, name(`msec', replace)
    twoway scatter mabsdiff pe, name(`absdiffc', replace)
    gr combine `msec' `absdiffc'
  }  
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





