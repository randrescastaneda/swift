/*==================================================
project:       
Author:        R.Andres Castaneda Aguilar 
----------------------------------------------------
Creation Date:    30 Jan 2020 - 18:39:10
==================================================*/

/*==================================================
0: Program set up
==================================================*/
program define swift_estimate, rclass
syntax varname [aw pw iw], ///
POVLine(numlist)           ///
Xvars(varlist)             ///
pe(real)                   /// pe in stepwise
pr(real)                   /// pr in stepwise
ADDm(integer)              /// number of m in mi 
seed(integer)              /// seed
bifold(varname)            ///
[                          ///
welfimp(string)            ///
povimp(string)            ///
]

version 16

marksample touse

//------------conditions
if ("`welfimp'" == "") {
  tempvar welfimp
}
if ("`povimp'" == "") {
  tempvar povimp
}

//------------Weigts
tempvar wi 
if missing(`"`weight'"') generate byte `wi' = 1
else generate `wi' `exp'
local wgtcall "[`weight' = `wi']"


//------------ right sample

local iff_tr "if (`touse' & `bifold' ==0)" // training

local iff_ts "if (`touse' & `bifold' ==1)" // test

//------------Main estimated

stepwise, pr(`pr') pe(`pe'): reg `varlist' `xvars' `wgtcall' `iff_tr'
local xvarm: colnames e(b)
local xvarm: subinstr local xvarm "_cons" "", word

// regression in training
reg `varlist' `xvarm' `wgtcall' `iff_tr'
return local r2_a = e(r2_a)

// Predict in test
tempvar yhat sqerr
predict `yhat' `iff_ts'
gen `sqerr' = (`varlist' - `yhat')^2 `iff_ts'
svy: mean `sqerr' `iff_ts'
return local mse = e(b)[1,1]

// set to missing all welfare values for mi
* replace `varlist' =. `iff_ts'

gen `welfimp' = `varlist' `iff_tr'
replace `welfimp' = . `iff_ts'


// Multiple imputation
mi set mlong                       // define output
mi register imputed `welfimp'      // register welfare variabe
mi impute reg `welfimp' `xvarm' `wgtcall', add(`addm') rseed(`seed') force

// new variabe
mi passive: gen `povimp' = (`welfimp' < `povline' )
mi estimate, post: svy: mean `povimp' `iff_ts'
return local mipov = e(b)[1,1]

end
exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
