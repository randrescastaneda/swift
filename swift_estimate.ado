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
bifold(varname)              

version 16

marksample touse

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
replace `varlist' =. `iff_ts'

// Multiple imputation
mi set mlong                       // define output
mi register imputed `varlist'      // register welfare variabe
mi impute reg `varlist' `xvarm' `wgtcall', add(`addm') rseed(`seed') force

// new variabe
tempvar imp_w
gen `imp_w' = (`varlist' < `povline' ) & _mi_m > 0
mi estimate, post: svy: mean `imp_w' `iff_ts'
return local mipov = e(b)[1,1]

end
exit
/* End of do-file */

><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
