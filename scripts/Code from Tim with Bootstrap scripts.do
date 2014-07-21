*===============================================================================
*
*
*    Across Data Set Wage Comparisons
*    (Runs matched regressions across many data sets for GED work)
*
*    By John Eric Humphries and Tim Kautz 8-2-2011
*
*
*
*===============================================================================


*===============================================================================
*
* Setup and Preparation
*
*==============================================================================
clear all
set mem 4g
set matsize 11000
pause on



local allclear = 1
local nlsy79   = 1  
local graphs   = 1



local ROOT_DIR "C:/Users/Tim/Documents/Research Projects/GED/Chapter 4/cross_section_graphs"
local FIG_DIR  "`ROOT_DIR'/figures"
local TAB_DIR  "`ROOT_DIR'/tables"
local INPUT_DIR  "`ROOT_DIR'/input"

local NLSY79data "`INPUT_DIR'/NLSY79 v8_5 - Partially Rebuilt Data to 2008 - Added Var Construction - Has Factors"


*=================================================================================
*
* Defining The Matices
*
*=================================================================================
*Creating Matrices

local matrix_rows=0
foreach var in LvlW LvlH LvlHrW Empl{
	foreach race in all {
		foreach sex in male female{
			foreach age in 2024 2529 3034 3539{
				foreach stub in AbilBg{
					foreach pse in AC {
						foreach edu in ged hsg{
							local matrix_rows=`matrix_rows'+1
						}
					}
				}
			}
		}
	}
}

local colnames "outcome race sex age control pse edu "

foreach var in ind_total dir_total tot_effect ind_nocoll_now b_nocoll_now p_nocoll_now ind_smcoll_now b_smcoll_now p_smcoll_now ind_degree_now b_degree_now p_degree_now {
	local colnames "`colnames' b_`var' se_`var' p_`var'"
}

local colnum : word count `colnames'

matrix define estimates=J(`matrix_rows',`colnum', .)
matrix colnames estimates = `colnames'





*=================================================================================
*
* NLSY79
*
*=================================================================================

use "`NLSY79data'" , clear

di "NSLY79"
di "NLSY79"

replace LvlHrW=. if LvlHrW==0
replace LvlH=. if LvlH==0

local outcomes "LvlW"
local background "mhgc urban_age14 faminc79 famstatus_broken  south_age14"
local standard_controls "dregion2 dregion3 dregion4 age black hisp  "
local standard_restrictions "  jail_yet!=1  & nowenrolled!=1  "  /* & (LvlH <=4000 | LvlH ==.) & (LvlW<300000 | LvlW==.) & LvlHrW>=3 & LvlHrW<=200 */
local ability " afqt_pre_std"
local hgc   "hgc_sec"
                                              
***********************************************************
* This runs all the cross-sectional regressions we graph
***********************************************************




* #3 Currently a GED, no college
gen		ged_nocoll_now = 1 if (dropout_now==0 & ged_now==1 & hsg_now==0 & touchcoll_now ==0)
replace ged_nocoll_now = 0 if ged_nocoll_now==.

* #4 Currently a high school grad, no college
gen		hsg_nocoll_now = 1 if (dropout_now==0 & ged_now==0 & hsg_now==1 & touchcoll_now ==0)
replace hsg_nocoll_now = 0 if hsg_nocoll_now==.

* #6 College dropout with a GED 
gen 	ged_smcoll_now = 1 if touchcoll_now==1 & aa_now==0 & gecoll_now==0 &college_nowenrolled!=1 & ged_now==1
replace ged_smcoll_now = 0 if ged_smcoll_now==.

* #7 College dropout with a HS Diploma
gen 	hsg_smcoll_now = 1 if touchcoll_now==1 & aa_now==0 & gecoll_now==0 &college_nowenrolled!=1 & hsg_now==1
replace hsg_smcoll_now = 0 if hsg_smcoll_now==.

*#8 GED with AA or BA
gen 	ged_degree_now = 1 if (aa_now==1 | gecoll_now==1) & ged_now==1
replace ged_degree_now = 0 if ged_degree==.

*#8 HSG with AA or BA
gen 	hsg_degree_now = 1 if (aa_now==1 | gecoll_now==1) & hsg_now==1
replace hsg_degree_now = 0 if hsg_degree==.


egen no_var=rowtotal(dropout_now ged_nocoll_now hsg_nocoll_now ged_smcoll_now hsg_smcoll_now ged_degree_now hsg_degree_now)
drop if no_var==0

local edu_regs dropout_now ged_nocoll_now hsg_nocoll_now ged_smcoll_now hsg_smcoll_now ged_degree_now hsg_degree_now

generate newid=id
tsset newid year




program define indirect_benefit, rclass
	
	syntax varlist [,if_statement(string)] 

	
	foreach var in ged_nocoll_now ged_smcoll_now ged_degree_now{
		mean `var' if `if_statement' & ged_now==1
		local p_`var'=_b[`var']
	}
	
	foreach var in hsg_nocoll_now hsg_smcoll_now hsg_degree_now{
		mean `var'  if `if_statement' & hsg_now==1
		local p_`var'=_b[`var'] 
	}
	
	reg `varlist' if `if_statement'
	
	foreach var in ged_nocoll_now ged_smcoll_now ged_degree_now hsg_nocoll_now hsg_smcoll_now hsg_degree_now{
		return scalar ind_`var'=_b[`var']*`p_`var''
		return scalar b_`var'=_b[`var']
		return scalar p_`var'=`p_`var''
	}
	
	return scalar ind_hsg_total=`p_hsg_smcoll_now'*_b[hsg_smcoll_now]+`p_hsg_degree_now'*_b[hsg_degree_now]
	return scalar dir_hsg_total=`p_hsg_nocoll_now'*_b[hsg_smcoll_now]
	return scalar tot_hsg_effect=`p_hsg_smcoll_now'*_b[hsg_smcoll_now]+`p_hsg_degree_now'*_b[hsg_degree_now]+`p_hsg_nocoll_now'*_b[hsg_smcoll_now]
	return scalar ind_ged_total=`p_ged_smcoll_now'*_b[ged_smcoll_now]+`p_ged_degree_now'*_b[ged_degree_now]
	return scalar dir_ged_total=`p_ged_nocoll_now'*_b[ged_nocoll_now]
	return scalar tot_ged_effect=`p_ged_smcoll_now'*_b[ged_smcoll_now]+`p_ged_degree_now'*_b[ged_degree_now]+`p_ged_nocoll_now'*_b[ged_nocoll_now]
end

local boot_return "ind_hsg_total=r(ind_hsg_total) dir_hsg_total=r(dir_hsg_total) tot_hsg_effect=r(tot_hsg_effect) ind_ged_total=r(ind_ged_total) dir_ged_total=r(dir_ged_total) tot_ged_effect=r(tot_ged_effect)"

foreach var in ged_nocoll_now ged_smcoll_now ged_degree_now hsg_nocoll_now hsg_smcoll_now hsg_degree_now{
	local boot_return "`boot_return' ind_`var'=r(ind_`var') b_`var'=r(b_`var') p_`var'=r(p_`var')"
}




local count=0

foreach outcome in `outcomes'  {
foreach race in all{
	if "`race'" == "all" local samplerest "& CrossSectSample==1"
	if "`race'" != "all" local samplerest ""
	
foreach pse in AC {	
foreach sex in  male  {
	
	local reg  "reg"

foreach age in 2024 2529 3034 3539 {
foreach control in Abil{
if "`pse'"=="AC" local ifstatement "  `standard_restrictions' & `sex'==1&A`age'==1 & `race'==1 `samplerest' "
	
	if "`control'"=="Raw"  		local controls `standard_controls'
	if "`control'"=="Abil"      local controls `ability'  `standard_controls'
	if "`control'"=="AbilBg"  	local controls `ability' `background' `standard_controls'
		
			
		
	if "`outcome'"=="LvlW" 			local outcome_num=1
	if "`outcome'"=="LvlH" 			local outcome_num=2
	if "`outcome'"=="Empl"	 		local outcome_num=3
	if "`outcome'"=="LvlHrW" 		local outcome_num=4
		
	if "`race'"=="all" 			local race_num=1
	if "`race'"=="white" 		local race_num=2
	if "`race'"=="black" 		local race_num=3
	if "`race'"=="hisp" 		local race_num=4
		
	if "`pse'"=="AC" 			local pse_num=1
	if "`pse'"=="NC" 			local pse_num=2
	if "`pse'"=="OC" 			local pse_num=3
		
	if "`sex'"=="male" 			local sex_num=1
	if "`sex'"=="female" 		local sex_num=2
		
	if "`age'"=="2024" 			local age_num=1
	if "`age'"=="2529" 			local age_num=2
	if "`age'"=="3034" 			local age_num=3
	if "`age'"=="3539" 			local age_num=4
	
	if "`control'"=="Raw" 		local control_num=1
	if "`control'"=="Abil" 		local control_num=2
	if "`control'"=="AbilBg" 	local control_num=3
	

	disp "`ifstatement'"
	bootstrap `boot_return', reps(10) seed(12345) cluster(id) idcluster(newid): indirect_benefit `outcome' ged_nocoll_now hsg_nocoll_now ged_smcoll_now hsg_smcoll_now ged_degree hsg_degree `controls', if_statement(`ifstatement')


	

	
	foreach hsedu in ged hsg{
		local count=`count'+1
		if "`hsedu'"=="ged" matrix estimates[`count',colnumb(estimates,"edu")]=1
		if "`hsedu'"=="hsg" matrix estimates[`count',colnumb(estimates,"edu")]=2
	
		foreach var in outcome race sex age control pse{
			matrix estimates[`count',colnumb(estimates,"`var'")] =``var'_num'
		}
		
		foreach est in ind_`hsedu'_total dir_`hsedu'_total tot_`hsedu'_effect ind_`hsedu'_nocoll_now b_`hsedu'_nocoll_now p_`hsedu'_nocoll_now ind_`hsedu'_smcoll_now b_`hsedu'_smcoll_now p_`hsedu'_smcoll_now ind_`hsedu'_degree_now b_`hsedu'_degree_now p_`hsedu'_degree_now ind_`hsedu'_nocoll_now b_`hsedu'_nocoll_now p_`hsedu'_nocoll_now ind_`hsedu'_smcoll_now b_`hsedu'_smcoll_now p_`hsedu'_smcoll_now ind_`hsedu'_degree_now b_`hsedu'_degree_now p_`hsedu'_degree_now{
			test _b[`est']==0
			local edu_pos_start=strpos("`est'","`hsedu'")-1
			disp "1"
			local edu_pos_end=strpos("`est'","`hsedu'")+4
			disp "1"
			local est_abrev1=substr("`est'",1,`edu_pos_start')
			disp "1"
			local est_abrev2=substr("`est'",`edu_pos_end',.)
			disp "1"
			local pos_string="`est_abrev1'`est_abrev2'"
			disp "`pos_string'"
	
			matrix estimates[`count',colnumb(estimates,"b_`pos_string'")]=_b[`est']
			matrix estimates[`count',colnumb(estimates,"se_`pos_string'")]=_se[`est']
			matrix estimates[`count',colnumb(estimates,"p_`pos_string'")]=r(p)
		}
	}
	
	matlist estimates
	
	
	
}


display "NLSY79 `var' `race' `sex' Age group `agegrp' "


} /* end agegrp loop */
} 
}
}
}



drop _all

svmat estimates

local count=0

foreach column in `colnames'{
	local count=`count'+1
	rename estimates`count' `column'
}

foreach column in outcome race sex age control pse{
	gen `column'_str = ""
}




	replace outcome_str="LvlW" 		if outcome==1
	replace	outcome_str="LvlH" 		if outcome==2
	replace outcome_str="Empl" 		if outcome==3
	replace outcome_str="LvlHrW" 	if outcome==4
	
	replace race_str  ="all" 		if race ==1
	replace race_str  ="white" 		if race ==2
	replace race_str  ="black" 		if race ==3
	replace race_str  ="hisp" 		if race ==4
		
	replace pse_str  ="AC" 			if pse ==1
	replace pse_str  ="NC" 			if pse ==2
	replace pse_str  ="OC" 			if pse ==3
		
	replace sex_str  ="male" 		if sex ==1
	replace sex_str  ="female" 		if sex ==2
		
	replace age_str  ="2024" 		if age ==1
	replace age_str  ="2529" 		if age ==2
	replace age_str  ="3034" 		if age ==3
	replace age_str  ="3539" 		if age ==4
	
	replace control_str  ="Raw" 	if control ==1
	replace control_str  ="Abil" 	if control ==2
	replace control_str ="AbilBg" 	if control ==3
	
foreach column in outcome race sex age control pse{
	drop `column'
	rename `column'_str `column'
}






*=================================================================================
*
*
*  THIS SECTION MAKES THE GRAPHS
*
*=================================================================================



drop if edu==.

gen num=_n

reshape long b_ se_ p_ , i(num ) j(state) string


rename b_ beta
rename  se_ se
rename p_ pvalue

gen model=.
local model=0

		
foreach age in 2024 2529{
	foreach edu in 1 2{
		foreach state in tot_effect ind_nocoll_now  ind_smcoll_now ind_degree_now {
			local model=`model'+1
			replace model=`model' if age=="`age'"&edu==`edu' & state=="`state'" 
		}
		local model=`model'+.5
	}
	local model=`model'+1
}

local model=0

foreach age in 3034 3539{
	foreach edu in 1 2{
		foreach state in tot_effect ind_nocoll_now  ind_smcoll_now ind_degree_now {
			local model=`model'+1
			replace model=`model' if age=="`age'"&edu==`edu' & state=="`state'" 
		}
		local model=`model'+.5
	}
	local model=`model'+1
}

gen upper=beta+se
gen lower=beta-se
	
gen 	sig1=beta if pvalue<=.05
replace sig1=. if pvalue>.05

gen 	sig2=beta if pvalue>.05&pvalue<=.1
replace sig2=. if pvalue<=.05|pvalue>.1

	
* gen min=.
* gen max=.
	
* foreach outcome in LvlH LvlHrW LvlW Empl{
		* foreach pse in all white black hisp{
			* foreach control in Raw Abil AbilBg{
				* egen max_temp=max(upper) if outcome=="`outcome'"&race=="`race'"&control=="`control'"
				* replace max=max_temp if outcome=="`outcome'"&race=="`race'"&control=="`control'"
				* drop max_temp
				* egen min_temp=min(lower) if outcome=="`outcome'"&race=="`race'"&control=="`control'"
				* replace min=min_temp if outcome=="`outcome'"&race=="`race'"&control=="`control'"
				* drop min_temp
		* }
	* }
* }

* replace min=0 if min>0

*  LvlH LvlHrW female 
gen beta_cur=.
gen se_cur=.
gen upper_cur=.
gen lower_cur=.
gen sig2_cur=.
gen sig1_cur=.

gen min_obs=.
gen max_obs=.

foreach outcome in LvlW {
	foreach sex in male {
		foreach  race in all {
			foreach control in Abil{
				replace min_obs=0
				replace max_obs=.
				
				egen max_obs_temp=max(upper) if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'"& (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")
				replace max_obs=max_obs_temp
				drop max_obs_temp
				
				egen min_obs_temp=min(lower) if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'"& (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")
				replace min_obs=min_obs_temp
				drop min_obs_temp
			
				foreach agegrp in A2029 A3039{
				if "`agegrp'"=="A2029"{
					replace beta_cur=beta if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
					
					replace upper_cur=upper if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
					
					replace lower_cur=lower if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
					
					replace se_cur=lower if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
					
					replace sig2_cur=sig2 if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
					
					replace sig1_cur=sig1 if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="2024"|age=="2529")
				}
				
				if "`agegrp'"=="A3039"{
					replace beta_cur=beta if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
					
					replace upper_cur=upper if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
					
					replace lower_cur=lower if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
					
					replace se_cur=lower if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
					
					replace sig2_cur=sig2 if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
					
					replace sig1_cur=sig1 if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" & (state=="tot_effect"|state=="ind_nocoll_now"|state=="ind_smcoll_now"|state=="ind_degree_now")&(age=="3034"|age=="3539")
				}
				
				
				twoway ///
				(bar beta_cur model if state=="tot_effect"&sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'",  color(white) lcolor(black) lwidth(.25)) ///
				(bar beta_cur model if state=="ind_nocoll_now"&sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'",  color(gs14) lcolor(black) lwidth(.25)) ///
				(bar beta_cur model if state=="ind_smcoll_now"&sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'",  color(gs6) lcolor(black) lwidth(.25)) ///
				(bar beta_cur model if state=="ind_degree_now"&sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'",  color(black) lcolor(black) lwidth(.25)) ///
				(rcap upper_cur lower_cur model if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'"  , color(gs0) lcolor(black) lwidth(.1)) ///
				(sc sig2_cur model if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'" , msymbol(o) mlcolor(black) mfcolor(white) msize(large) ) ///
				(sc sig1_cur model if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'"  , msymbol(o) mcolor(black) msize(large) ) ///
				(sc min_obs model  if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'", msymbol(o) mcolor(none) msize(tiny) ) ///
				(sc max_obs model  if sex=="`sex'"&outcome=="`outcome'"&race=="`race'"&control=="`control'", msymbol(o) mcolor(none) msize(tiny) ) ///
				, scheme(s1manual) legend( rows(2) order( 1 "Total" 2 "Direct" 7 "p<0.05(vs.0)" 5  "S.E."  3  "Some Col." 4 "AA/BA"  6 "p<0.10(vs.0)"  )) yline(0, lcolor(black) ) xtitle("", size(small)) `ylabel' `yscale' title("")  xlabel( 4 "Age 20 to 24" 15 "Age 25 to 29" ,  angle(0)  labsize(medlarge) noticks custom labgap(*8) )  ///
				xlabel( 2.5 "GED" 7 "HSG" 12.5 "GED" 17 "HSG" ,  angle(0)  labsize(medsmall) noticks custom add )
				local ROOT_DIR "C:/Users/Tim/Documents/Research Projects/GED/Chapter 4/cross_section_graphs"
				local FIG_DIR  "`ROOT_DIR'/figures"
				graph export "`FIG_DIR'/PSE_NLSY79ByAge_`outcome'_`sex'_`race'_`control'_`agegrp'.eps" , replace logo(off)
				graph export "`FIG_DIR'/PSE_NLSY79ByAge_`outcome'_`sex'_`race'_`control'_`agegrp'.png" , replace 
				pause
				
				foreach stat in beta_cur sig2_cur sig1_cur se upper_cur lower_cur{
					replace `stat'=.
				}
				
			}
		}
	}
}
}	 







