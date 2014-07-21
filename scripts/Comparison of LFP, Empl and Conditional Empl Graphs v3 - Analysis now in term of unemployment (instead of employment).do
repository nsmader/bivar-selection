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
set more off, perm
pause on

local RunEstimation = 0
local RunGraphs     = 1


local ROOT_DIR "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Cross-Sectional Estimates/"
local OUTPUT_DIR  "`ROOT_DIR'/output"
local FIG_DIR  "`ROOT_DIR'/output/figures"
local TAB_DIR  "`ROOT_DIR'/output/tables"

local NLSY79data "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Intermediate Data/NLSY79 v8_5 - Partially Rebuilt Data to 2008 - Added Var Construction - Has Factors"


*=================================================================================
*
* Defining The Matices
*
*=================================================================================
*Creating Matrices

local matrix_rows=0
foreach var in LFP Empl UnempIfLFP {
	foreach race in all white black hisp {
		foreach sex in male female {
			foreach age in 2024 2529 3034 3539 {
				foreach stub in Raw Abil AbilBg {
					foreach pse in AC NC { /* OC */
						foreach edu in ged hsg {
							local matrix_rows=`matrix_rows'+1
						}
					}
				}
			}
		}
	}
}

local ColNames = "outcome race sex age control pse edu beta se sig1 sig2"
matrix define estimates=J(`matrix_rows', wordcount("`ColNames'"), .)
matrix colnames estimates = `ColNames'



*=================================================================================
*
* NLSY79
*
*=================================================================================

if `RunEstimation' == 1 {

	use "`NLSY79data'" , clear

	replace LvlHrW = . if LvlHrW == 0
	replace LvlH = . if LvlH == 0
	replace Empl = 0 if LFP == 0
	gen EmplIfLFP = Empl if LFP == 1
	gen UnempIfLFP = Unempl if LFP == 1

	local outcomes              " LFP Empl UnempIfLFP"
	local background            " mhgc faminc79 famstatus_broken urban_age14 south_age14"
	local standard_controls     " dregion2 dregion3 dregion4 _Iyear* age black hisp  "
	local standard_restrictions " & jail_yet!=1  & enrolled!=1  "  /* & (LvlH <=4000 | LvlH ==.) & (LvlW<300000 | LvlW==.) & LvlHrW>=3 & LvlHrW<=200 */
	local ability               " afqt_pre_std "
	local hgc   "" /* hgc_sec */
												  
	***********************************************************
	* This runs all the cross-sectional regressions we graph
	***********************************************************

	local count = 0

	foreach outcome in `outcomes'  {
	foreach race in all black white hisp {
		if "`race'" == "all" local samplerest "& CrossSectSample == 1"
		if "`race'" != "all" local samplerest ""
		
	foreach pse in AC NC {	/* OC */
	foreach sex in  male female  {
		if "`pse'"=="AC" local ifstatement "  if (ged_now==1 | hsg_now==1 | dropout_only==1)   `samplerest' `standard_restrictions' & `sex'==1 & `race'==1  "
		if "`pse'"=="NC" local ifstatement "  if (ged_only==1 | hsg_only==1 | dropout_only==1) `samplerest' `standard_restrictions' & `sex'==1 & `race'==1   "
		if "`pse'"=="OC" local ifstatement "  if ( (ged_now==1 & ged_only!=1) | (hsg_now==1 & hsg_only!=1) | dropout_only==1) `samplerest' `standard_restrictions' & `sex'==1 & `race'==1   "
		local reg  "reg"

	foreach age in 2024 2529 3034 3539 {
	foreach control in Raw Abil AbilBg {

		if "`control'"=="Raw"  		local controls `standard_controls'
		if "`control'"=="Abil"      local controls `ability' `standard_controls'
		if "`control'"=="AbilBg"  	local controls `ability' `background' `standard_controls'
			
		if "`outcome'"=="LFP" 		  local outcome_num=1
		if "`outcome'"=="Empl" 		  local outcome_num=2
		if "`outcome'"=="UnempIfLFP"  local outcome_num=3
			
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
			
		qui `reg' `outcome'  ged_now hsg_now  `controls' `ifstatement'  & A`age'==1, vce(cluster id)
		
		local sig_ged_drop=0
		quietly test (_b[ged_now]=0)
		local sig=0
		if r(p)<.05 local sig_ged_drop=1
		
		local sig_hsg_drop=0
		quietly test (_b[hsg_now]=0)
		local sig=0
		if r(p)<.05 local sig_hsg_drop=1
		
		local sig_ged_hsg=0
		quietly test (_b[ged_now]=_b[hsg_now])
		local sig=0
		if r(p)<.05 local sig_ged_hsg=1
		
		local count=`count'+1
		
		/* Save the (numerical) categorical value for value of each loop that's run... */
			/* ... for GEDs */
			foreach var in outcome race sex age control pse {
				matrix estimates[`count',colnumb(estimates,"`var'")] = ``var'_num'
			}
			
			matrix estimates[`count',colnumb(estimates,"edu")]  = 1
			matrix estimates[`count',colnumb(estimates,"beta")] = _b[ged_now]
			matrix estimates[`count',colnumb(estimates,"se")]   = _se[ged_now]
			
			matrix estimates[`count',colnumb(estimates,"sig1")] = `sig_ged_drop'
			matrix estimates[`count',colnumb(estimates,"sig2")] = `sig_ged_hsg'

			local count=`count'+1
			
			/* ... for HSGs */
			foreach var in outcome race sex age control pse {
				matrix estimates[`count',colnumb(estimates,"`var'")] =``var'_num'
			}
			
			matrix estimates[`count',colnumb(estimates,"edu")]  = 2
			matrix estimates[`count',colnumb(estimates,"beta")] = _b[hsg_now]
			matrix estimates[`count',colnumb(estimates,"se")]   = _se[hsg_now]
			matrix estimates[`count',colnumb(estimates,"sig1")] = `sig_hsg_drop'
			matrix estimates[`count',colnumb(estimates,"sig2")] = `sig_ged_hsg'
		
		display "NLSY79 outcome = `outcome', race = `race', pse = `pse', sex = `sex', age = `age', control = `control'"

	} /* control */
	} /* end agegrp loop */
	} /* sex */
	} /* pse */
	} /* race */
	} /* outcome */


	/*------------------------*/
	/* * * SAVE ESTIMATES * * */
	/*------------------------*/

		drop _all

		svmat estimates

		/* NSM: This saves the columns of the "estimates" matrix as new data (everything else was dropped by the "drop _all").
			Each column is generically called "estimates#", so the following code renames it. */

		local count = 0

		foreach column in outcome race sex age control pse edu beta se sig1 sig2 {
			local count = `count'+1
			rename estimates`count' `column'
		}

		/* Note: we leave it as acceptable that edu \in {1,2}, and the beta, se, sig1 and sig2 values are as they were when they were saved */

		foreach column in outcome race sex age control pse {
			gen `column'_str = ""
		}
			replace outcome_str = "LFP" 	    if outcome==1
			replace	outcome_str = "Empl" 	    if outcome==2
			replace outcome_str = "UnempIfLFP" if outcome==3
			
			replace race_str = "all" 		  if race ==1
			replace race_str = "white" 		  if race ==2
			replace race_str = "black" 		  if race ==3
			replace race_str = "hisp" 		  if race ==4
				
			replace pse_str  = "AC" 		  if pse ==1
			replace pse_str  = "NC" 		  if pse ==2
			replace pse_str  = "OC" 		  if pse ==3
				
			replace sex_str  = "male" 		  if sex ==1
			replace sex_str  = "female" 	  if sex ==2
				
			replace age_str  = "2024" 		  if age ==1
			replace age_str  = "2529" 		  if age ==2
			replace age_str  = "3034" 		  if age ==3
			replace age_str  = "3539" 		  if age ==4
			
			replace control_str = "Raw" 	  if control ==1
			replace control_str = "Abil" 	  if control ==2
			replace control_str = "AbilBg" 	  if control ==3
			
		foreach column in outcome race sex age control pse {
			drop `column'
			rename `column'_str `column'
		}

	save "`OUTPUT_DIR'/Full_Reg_Estimates.dta", replace

} /* Finished Running Estimation */

*=================================================================================
*
*
*  THIS SECTION MAKES THE GRAPHS
*
*=================================================================================

if `RunGraphs' == 1 {

	use "`OUTPUT_DIR'/Full_Reg_Estimates.dta", clear

	local count = 0
	gen model = .
	foreach age in 2024 2529 3034 3539 {
		/* NSM: Want to redefine this loop to be my inner display piece in the graph. For me, this is across outcomes rather than races. */
		foreach outcome in LFP Empl UnempIfLFP {
			local count = `count' + 1
			replace model = `count' if age=="`age'" & outcome=="`outcome'" & edu==1
			local count = `count' + 1
			replace model = `count' if age=="`age'" & outcome=="`outcome'" & edu==2
		}
		local count = `count' + 1
	}

		gen upper = beta + se
		gen lower = beta - se
		
		replace sig2 = beta*sig2
		replace sig1 = beta*sig1
		replace sig2 = . if sig2==0
		replace sig1 = . if sig1==0

		gen min=.
		gen max=.
		
	foreach outcome in LFP Empl UnempIfLFP {
			foreach pse in AC NC { /* OC */
				foreach control in Raw Abil AbilBg {
					egen max_temp = max(upper) if outcome == "`outcome'" & pse == "`pse'" & control == "`control'"
					replace max = max_temp     if outcome == "`outcome'" & pse == "`pse'" & control == "`control'"
					drop max_temp
					egen min_temp = min(lower) if outcome == "`outcome'" & pse == "`pse'" & control == "`control'"
					replace min = min_temp     if outcome == "`outcome'" & pse == "`pse'" & control == "`control'"
					drop min_temp
			}
		}
	}

	foreach race in all white black hisp {
		foreach sex in male female {
			foreach pse in AC NC { /* OC */
				foreach control in Raw Abil AbilBg {
					#delimit;
					twoway (bar beta model if edu==1 & sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'", color(gs4)  lcolor(black) lwidth(.25))
						   (bar beta model if edu==2 & sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'", color(gs14) lcolor(black) lwidth(.25))
					(rcap upper lower model if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'", color(gs0) lcolor(black) lwidth(.1))
					(sc sig2 model if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'" & edu==1, msymbol(o) mlcolor(black) mfcolor(white) msize(huge))
					(sc sig1 model if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'" & edu==1, msymbol(o) mcolor(black) msize(large))
					(sc sig1 model if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'" & edu==2, msymbol(d) mcolor(black) msize(large))
					(sc min model  if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'", msymbol(o) mcolor(none) msize(tiny) )
					(sc max model  if sex=="`sex'" & race=="`race'" & pse=="`pse'" & control=="`control'", msymbol(o) mcolor(none) msize(tiny) ),
					legend( rows(2) order( 1 "GED" 4 "p<0.05 (GED vs.HSG)" 5 "p<0.05 (GED vs.Drop)" 2 "HSG" 6 "p<0.05 (HSG vs.Drop)" 3 "S.E."  ) symxsize(*.75)) yline(0, lcolor(black) ) xtitle("", size(small))
					xlabel( 3.5 "Age 20 to 24" 10.5 "Age 25 to 29"  17.5 "Age 30 to 34"  24.5 "Age 35 to 39"  ,  angle(0)  labsize(medlarge) noticks custom labgap(*15) ) xlabel( 1.5 "LFP" 3.5 "Empl" 5.5 "Unemp" 8.5 "LFP" 10.5 "Empl" 12.5 "Unemp" 15.5 "LFP" 17.5 "Empl" 19.5 "Unemp" 22.5 "LFP" 24.5 "Empl" 26.5 "Unemp",  angle(-45)  labsize(medsmall) noticks custom add )
					scheme(s1manual)  `ylabel' `yscale' title("");
					graph export "`FIG_DIR'/LFPReg_NLSY79ByAge_`race'_`sex'_`pse'_`control'.eps" , replace logo(off) fontface(Times);
					graph export "`FIG_DIR'/LFPReg_NLSY79ByAge_`race'_`sex'_`pse'_`control'.png" , replace;
					#delimit cr;
				} /* control */
			} /* pse */
		} /* sex */
	} /* race */
	 
} /* Finished Running Graphs */


