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

qui foreach i of numlist 1/`nsplit' {
  
  // sample conditions
  local iff_tr "if (`touse' & `fold' !=`i')" // training
  local iff_ts "if (`touse' & `fold' ==`i')" // test
  
  // Actual poverty rate in fold i
  svy: mean `pov' `iff_ts'
	local rpov = e(b)[1,1]	
  
  
  local pe = `stpe'
  local pr = `pe' + `tolerancepe'
  
  while `pe' <= `limitpe' {
		use `basef', clear
		
		// stepwise, get variables
		stepwise, pr(`pr') pe(`pe'): reg `varlist' `xvars' `wgtcall' `iff_tr'
		local xvars: colnames e(b)
		local xvars: subinstr local xvars "_cons" "", word
		
		// regression in training
		reg `varlist' `xvars' `wgtcall' `iff_tr'
    local r2_a = e(r2_a)[1,1]
		
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
		mi impute reg `varlist' `xvars' `wgtcall', add(`addm') rseed(`seed') force
		
		// new variabe
		tempvar imp_w
		gen `imp_w' = (`varlist' < `povline' ) & _mi_m > 0
		mi estimate, post: svy: mean `imp_w' `iff_ts'
    local mipov = e(b)[1,1]
    
    local diffmse = abs(`rpov' `mipov')
		
		matrix `X' = nullmat(`X') \ `pe', `i', `rpov', `r2_a', `mse', `mipov', `diffmse'
    
    
    local pe = `pe' + `deltape'
    local pr = `pe' + `tolerancepe'
  }
	
}

/*==================================================
3: Organize results
==================================================*/

//------------return matrix X
mat colnames `X' = pe fold_ poor_ r2 mse imp_poor absdiff
return matrix pmi = `X'

//------------calculate mean of MSE  and diff in poverty rates by pe
drop _all 
svmat `X', n(col)
collapse (mean) mmse = mse (mean) mabsdiff = absdiff, by(pe)

tempname Y
mkmat *, matrix(`Y')

return matrix mmse = `Y'

if ("`chart'" != "nochart") {
  tempfile msec absdiffc 
  twoway scatter mmse pe
  gr save `msec', replace
  
  twoway scatter mabsdiff pe
  gr save `absdiffc', replace
  
  gr combine `msec' `absdiffc'
}

end
exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

Notes:
1.
2.
3.


Version Control:


