/*=================================================================================================================*/
/*                                                                                                                 */
/*    Generate Graphs Comparing Estimation for NLSY79 Including Selection Correction for Labor Force Participation */
/*                                                                                                                 */
/*    By Nicholas Mader 12/8/2011                                                                                  */
/*                                                                                                                 */
/*=================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------*/


/*======================================*/
/*                                      */
/* * * SET ENVIRONMENTAL PARAMETERS * * */
/*                                      */
/*======================================*/

	clear all
	cap set mem 8g
	cap set maxvar 30000
	cap set matsize 10000
	set more off, perm
	pause on
	set trace off
	set tracedepth 1
	
	local RunEstimation = 0
	local RunGraphing   = 1

	local INPUT_DIR  "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Intermediate Data/"
	local ROOT_DIR   "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Cross-Sectional Estimates"
	local FIG_DIR    "`ROOT_DIR'/output/figures"
	local TAB_DIR    "`ROOT_DIR'/output/tables"
	local MAT_DIR    "`ROOT_DIR'/output/matrices"
	

	local NLSY79data "`INPUT_DIR'/NLSY79 v8_5 - Partially Rebuilt Data to 2008 - Added Var Construction - Has Factors.dta"
	/*
		use "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Intermediate Data//NLSY79 v8_5 - Partially Rebuilt Data to 2008 - Added Var Construction - Has Factors.dta", clear
	*/

	/*-----------------------------------------------------------------------*/
	/* Define Program for Bootstrapping Parameters in the Selection Equation */
	/*-----------------------------------------------------------------------*/
	
		cap program drop SelectionEst
		program define SelectionEst, rclass
		
			syntax varlist, control_vars(string) sel_vars(string) if_statement(string)
		
			/* Attempt 1 (successful) 
				sum  afqt_pre_std if LFP == 1
				return scalar mymean = r(mean) */
				
			/* Attempt 2 (successful) 
				reg `varlist' `control_vars' `if_statement' & LFP == 1, vce(cluster id)
				return scalar b_ged_now = _b[ged_now]
				return scalar b_hsg_now = _b[hsg_now]*/


			/*	Attempt 3 (pending) */
			probit LFP `control_vars' `sel_vars' `if_statement'
				cap drop zg zg_p zg_c invmills
				predict zg, xb
				gen zg_p = normalden(zg)
				gen zg_c = normal(zg)
				gen invmills = -(zg_p/zg_c)
			
			di "reg `varlist' `control_vars' invmills `if_statement'"
			preserve
			keep if LFP == 1
			reg `varlist' `control_vars' invmills `if_statement'
			return scalar b_ged_now = _b[ged_now]
			return scalar b_hsg_now = _b[hsg_now]
			restore
			
		end



/*============================*/
/*                            */
/* * * SET RUN PARAMETERS * * */
/*                            */
/*============================*/
	
	local DepVarList    = " LvlW Empl LvlHrW LvlH" /* Desc */
	local RaceVarList   = "all white black hisp"
	local GenderVarList = "male female"
	local PseRunList    = "NC AC"
	local AgeGroupList  = "2024 2529 3034 3539"
		local nAge = wordcount("`AgeGroupList'")
		local Age2024Label = "Ages 20-24"
		local Age2529Label = "Ages 25-29"
		local Age3034Label = "Ages 30-34"
		local Age3539Label = "Ages 35-39"
		
	local LFPSpecList   = "Raw Abil AbilBg" /*  */
		local  nStage1Spec = wordcount("`LFPSpecList'")
	local SpecEstList   = "AbilBg Cond Selec" /* Raw */
		local  nRegSpec    = wordcount("`SpecEstList'")
	local SpecGraphList = "AbilBg Cond Selec" /* Raw  */
		local RawLabel    = "Raw"
		local AbilLabel   = "Abil"
		local AbilBgLabel = "AbBg"
		local CondLabel   = "Cond"
		local SelecLabel  = "Sel"

	local BasicControls   = "dregion2 dregion3 dregion4 _Iyear* black hisp" /* jail_yet TimeInJail */
	local BgControls      = "mhgc urban_age14 faminc79 famstatus_broken south_age14 "
	local AbilityControls = "afqt_pre_std"
	local hgc             = " " /* "hgc_sec_ever" */
	local SelVars         = "unemp1 unemp2 unemp3 unemp4  lowage1 lowage2 lowage3 lowage4  married_now num_children_in_hh baby_in_hh toddler_in_hh SpouseWage "
	local SampleRestr     = "jail_yet!=1 & enrolled!=1 "  /* & (LvlH <=4000 | LvlH ==.) & (LvlW<300000 | LvlW==.) & LvlHrW>=3 & LvlHrW<=200 */
	

/*========================*/
/*                        */
/* * * RUN ESTIMATION * * */
/*                        */
/*========================*/

if `RunEstimation' == 1 {

	use "`NLSY79data'" , clear
	cap gen newid = id
	  
	keep `DepVarList' `BgControls' LFP `BasicControls' `AbilityControls' `SelVars' `hgc' ///
		id male female white all ged_now ged_only dropout_now dropout_only hsg_only hsg_now jail_now jail_yet enrolled CrossSectSample A2024 A2529 A3034 A3539
		
	/* Make small data adjustments */
	
		replace Empl   = 0 if LFP == 0
		replace LvlHrW = . if LFP == 0
		replace LvlH   = . if LFP == 0
		
		gen     EdCat = 1 if dropout_now == 1
		replace EdCat = 2 if ged_now == 1
		replace EdCat = 3 if hsg_now == 1
		label define EdCatLabel 1 "d" 2 "g" 3 "h"
		label values EdCat EdCatLabel
	

	foreach DepVar in `DepVarList' { 
	foreach GenderVar in `GenderVarList' {
		local Gx = "G" + substr("`GenderVar'",1,1)	
	foreach RaceVar in `RaceVarList' {
		local Rx = "R" + substr("`RaceVar'",1,1)
		if "`RaceVar'" == "all" local RaceCond "& CrossSectSample == 1"
		if "`RaceVar'" != "all" local RaceCond "& `RaceVar' == 1"
		
		eststo clear
		
	foreach Pse in `PseRunList' {
		if "`Pse'" == "AC" local EdCond " & (ged_now==1 | hsg_now==1 | dropout_now==1)    "
		if "`Pse'" == "NC" local EdCond " & (ged_only==1 | hsg_only==1 | dropout_only==1) "
		if "`Pse'" == "OC" local EdCond " & ( (ged_now==1 & ged_only!=1) | (hsg_now==1&hsg_only!=1) | dropout_only==1) "
	
	
	/*----------------------------------------*/
	/* Define the Matrices for Saving Results */
	/*----------------------------------------*/
		
		local  desc_mat      = "LFPDesc_`Gx'`Rx'`Pse'"
		local  nDescCols     = 4 /* This will be DO LFP, GED LFP, DO Empl, GED Empl */
		local  DescColNames  = "beta se sig dropout ged"
		matrix `desc_mat'    = J((`nDescCols' + 1)*`nAge' - 1, wordcount("`DescColNames'"), .)  /* We add an extra slot for a missing value which will separate different series of statistics by age */
		matrix colnames `desc_mat' = `DescColNames'
		
		local  Stage1_mat     = "LFPStage1_`Gx'`Rx'`Pse'"
		local  Stage1ColNames = "beta se sig ged hsg `LFPSpecList'"
		matrix `Stage1_mat'   = J((`nStage1Spec' + 1)*`nAge' - 1, wordcount("`Stage1ColNames'"), .)
		matrix colnames `Stage1_mat' = `Stage1ColNames'
		
		local  reg_mat     = "LFPReg_`DepVar'_`Gx'`Rx'`Pse'"
		local  RegColNames = "beta se sig1 sig2 ged hsg `SpecEstList'"
		matrix `reg_mat'   = J((`nRegSpec'*2+1)*`nAge' - 1, wordcount("`RegColNames'"), .)
		matrix colnames `reg_mat' = `RegColNames'
		
		local ixDesc   = 1
		local ixStage1 = 1
		local ixReg    = 1
	
	/*------------------------------------*/
	/* Run Loops for the Various Analyses */
	/*------------------------------------*/
	
	foreach AgeGrp in `AgeGroupList' {
	
		local ifstatement  = "if `SampleRestr' & A`AgeGrp' == 1 & `GenderVar'==1 `RaceCond' `EdCond' `SpecCond'"
		
			
		if "`DepVar'" == "Desc" {
		
		/*--------------------------------------------------------------------------------*/
		/* * * Collect Summary Statistics on LFP and Empl for Dropouts, GEDs and HSGs * * */
		/*--------------------------------------------------------------------------------*/
		
			#delimit;
			mean LFP Empl `ifstatement', over(EdCat);
			
			foreach v in LFP Empl {;
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"dropout")] = 1;
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"beta")]    = [`v']_b[d];
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"se")]      = [`v']_se[d];
					local ixDesc = `ixDesc' + 1;
				
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"ged")]     = 1;
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"beta")]    = [`v']_b[g];
				matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"se")]      = [`v']_se[g];
				test [`v']_b[d] = [`v']_b[g];
				if r(p) < 0.05 matrix `desc_mat'[`ixDesc', colnumb(`desc_mat',"sig")] = [`v']_b[g];
					local ixDesc = `ixDesc' + 1;
				
			};
			local ixDesc = `ixDesc' + 1; /* This adds an extra space inbetween blocks of statistics by age. */
	

		/*------------------------------------------------------------*/
		/* * * Collect Results from First Stage of LFP Estimation * * */
		/*------------------------------------------------------------*/

			#delimit;
			foreach Stage1Spec in `LFPSpecList' {;
				if "`Stage1Spec'" == "Raw"    local controls `BasicControls';
				if "`Stage1Spec'" == "Abil"   local controls `BasicControls' `AbilityControls';
				if "`Stage1Spec'" == "AbilBg" local controls `BasicControls' `BgControls' `AbilityControls';				
			
				eststo LFP_`Stage1Spec'_`Gx'`Rx'`Pse': qui probit LFP ged_now hsg_now `controls' `ifstatement', vce(cluster id);
				
				foreach v in ged hsg {;
				
					matrix `Stage1_mat'[`ixStage1', colnumb(`Stage1_mat',"`v'")]          = 1;
					matrix `Stage1_mat'[`ixStage1', colnumb(`Stage1_mat',"`Stage1Spec'")] = 1;
					matrix `Stage1_mat'[`ixStage1', colnumb(`Stage1_mat',"beta")]         = _b[`v'_now];
					matrix `Stage1_mat'[`ixStage1', colnumb(`Stage1_mat',"se")]           = _se[`v'_now];
					test _b[`v'_now] = 0;
						if r(p) < 0.05 matrix `Stage1_mat'[`ixDesc', colnumb(`Stage1_mat',"sig1")] = _b[`v'_now];
					test _b[ged_now] = _b[hsg_now];
						if r(p) < 0.05 matrix `Stage1_mat'[`ixDesc', colnumb(`Stage1_mat',"sig2")] = _b[`v'_now];
					local ixStage1 = `ixStage1' + 1;
				};
				
			};
			local ixStage1 = `ixStage1' + 1; /* This adds an extra space inbetween blocks of statistics by age. */
			#delimit cr;
			
		}
		
		
		/*-----------------------------*/
		/* * * Run Regression Code * * */
		/*-----------------------------*/
		
		if "`DepVar'" != "Desc" {
		
			foreach Spec in `SpecEstList' {
				if "`Spec'" == "Raw"    local controls `BasicControls'
				if "`Spec'" == "Abil"   local controls `BasicControls' `AbilityControls'
				if "`Spec'" == "AbilBg" local controls `BasicControls' `BgControls' `AbilityControls'
				if "`Spec'" == "Cond"   local controls `BasicControls' `BgControls' `AbilityControls'
				if "`Spec'" == "Selec"  local controls `BasicControls' `BgControls' `AbilityControls'
				
				if "`Spec'" == "Cond" local SpecCond = "& LFP == 1"
				if "`Spec'" != "Cond" local SpecCond = ""
				
				di "Running DepVar = `DepVar', GenderVar = `GenderVar', RaceVar = `RaceVar', Pse = `Pse', Spec = `Spec'"
		
			/*----------------------------------------*/
			/* * * Run Estimaton and Save Results * * */
			/*----------------------------------------*/
			
				#delimit;
				if "`Spec'" != "Selec" {;
					eststo R_`Spec'_`DepVar'_`AgeGrp'`Gx'`Rx'`Pse': qui reg `DepVar' ged_now hsg_now `controls' `ifstatement' `SpecCond', `reg_opts' vce(cluster id);
				};

				if "`Spec'" == "Selec" {; 
						bootstrap ged_est=r(b_ged_now) hsg_est=r(b_hsg_now), reps(100) seed(12345) idcluster(newid) cluster(id):
							SelectionEst `DepVar', control_vars(ged_now hsg_now `controls') sel_vars(`SelVars') if_statement(`ifstatement');
				};
				
				
					
				if "`Spec'" != "Selec" {; local b_ged  _b[ged_now]; local se_ged  _se[ged_now];
										  local b_hsg  _b[hsg_now]; local se_hsg  _se[hsg_now]; };
				if "`Spec'" == "Selec" {; local b_ged  _b[ged_est]; local se_ged  _se[ged_est];
										  local b_hsg  _b[hsg_est]; local se_hsg  _se[hsg_est]; };
				
				foreach v in ged hsg {;
				
					matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"`v'")]    = 1;
					matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"`Spec'")] = 1;
					matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"beta")]   = `b_`v'';
					matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"se")]     = `se_`v'';
					quietly test (`b_`v'' = 0);
						if r(p) < 0.05 matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"sig1")] = `b_`v'';
					quietly test (`b_ged' = `b_hsg');
						if r(p) < 0.05 matrix `reg_mat'[`ixReg', colnumb(`reg_mat',"sig2")] = `b_`v'';
					local ixReg = `ixReg' + 1;
					
				};

				#delimit cr;

			} /* End of Loop Across Regression Specifications */
			
			local ixReg = `ixReg' + 1
			
		} /* End of Estimation Code */

		display "Finished Running: DepVar = `DepVar', RaceVar = `RaceVar', GenderVar = `GenderVar',  AgeGrp = `AgeGrp', Pse = `Pse'."

	} /* End of AgeGrp Loop */
	
		matsave `reg_mat',    path("`MAT_DIR'") sav replace
		matsave `desc_mat',   path("`MAT_DIR'") sav replace
		matsave `Stage1_mat', path("`MAT_DIR'") sav replace
	
	} /* End of Pse Loop */
		esttab R_* using "`TAB_DIR'/NLSY79 Cross-Sectional Estimated Returns - Sensitivity to Control for LFP - `DepVar' - `RaceVar' `GenderVar's.csv", compress nostar replace csv nonumbers mtitles se
		if strpos("`DepVar'","Desc") >0 esttab LFP_* using "`TAB_DIR'/NLSY79 Cross-Sectional Estimated Returns - LFP First Stage - `RaceVar' `GenderVar's.csv", compress nostar replace csv nonumbers mtitles se
	} /* End of Race Loop */
	} /* End of Gender Loop */
	} /* End DepVar Loop */

} /* End of Estimation */

clear


/*=========================================================*/
/*                                                         */
/* * *  GGGGGRRRRRAAAAAPPPPHHHH MMMMAAAAKKKKEEERR!!!!! * * */
/*                                                         */
/*=========================================================*/

if `RunGraphing' == 1 {
	
	foreach Pse     in `PseRunList'  {
	foreach RaceVar in `RaceVarList' {
		local Rx = "R" + substr("`RaceVar'",1,1)
	foreach GenderVar in `GenderVarList' {
		local Gx = "G" + substr("`GenderVar'",1,1)
		
	/*-----------------------------*/
	/* Generate Descriptive Graphs */
	/*-----------------------------*/

	/*-------------------------------------*/
	/* Generate Graphs of First Stage Info */
	/*-------------------------------------*/
			
			
	/*----------------------------*/
	/* Graph Regression Estimates */
	/*----------------------------*/
		
	foreach DepVar  in `DepVarList'  {
		 
		/* Identify the Minimum and Maximum Estimates Across Age and Specification (but within Gender, given how disparate the estimates are) */
		
			local max_y = 0
			local min_y = 0
			
			foreach AgeGrp in `AgeGroupList' {
				foreach Spec in `SpecGraphList' {
					foreach edu in ged hsg {
						local reg_mat LFPReg_`DepVar'_`Gx'`Rx'`Pse'
						use "`MAT_DIR'/`reg_mat'.dta", clear
						gen upper = beta + se
						gen lower = beta - se
						qui sum upper
						local max_y = max(`max_y',r(max))
						qui sum lower
						local min_y = min(`min_y',r(min))
					} /* End of Loop Across Education */
				} /* End of Loop Across specifications */
			} /* End of Loop Across Ages */

		/*-----------------------------------------------*/
		/* Contruct Components of the Graphing Statement */
		/*-----------------------------------------------*/
			di "123"
			#delimit;
			use "`MAT_DIR'/`reg_mat'.dta", clear;
				* replace sig1 = . if sig1 == -999;
				* replace sig2 = . if sig2 == -999;
				/* *** Need to make sure that we select only the rows of analysis that we want *** */
				gen x = _n;
				gen upper = beta + se;
				gen lower = beta - se;
				gen min=`min_y';
				gen max=`max_y';
			
			/* Series Plot Statement and Legend Order Statement */
				local PlotStatement; local LegendOrder; local i = 1;
				local fcolor_ged gs4; local fcolor_hsg gs14;
				foreach e in ged hsg {;
					local PlotStatement `PlotStatement' (bar beta x if `e' == 1, lcolor(black) fcolor(`fcolor_`e'') lwidth(0.25)); 
					local LegendOrder `"`LegendOrder' `i' "``e'Label'" "'; /* " */ 
					local i = `i' + 1; };
				local ip1 = `i' + 1; local ip2 = `i' + 2;
				local LegendOrder `" `LegendOrder' `i' "S.E." `ip2' "p<0.05 (vs. DO)" `ip1' "p<0.05 (GED vs. HSG)" "';  /* " */
					/* Note that the ip# numbers are reversed because the GEDvsHSG dot is plotted first in order to be a background to the vsDO dot */
				di `" `PlotStatement' "'; /* " */
				di `" `Legend Order' "'; /* " */				
			
			/* X Label for Ages */
				local AgeXLabel;
				local a = (`nRegSpec'*2+1)/2; /* This centers the age label to go with groups of estimates */
				foreach AgeGroup in `AgeGroupList' {;
					local AgeXLabel `"`AgeXLabel' `a' "`Age`AgeGroup'Label'" "'; /* " */
					local a = `a' + (`nRegSpec'*2) + 1; };
				 di `"`AgeXLabel'"'; /* " */
			
			/* X Label for Specifications */
				local SpecXLabel;
				local s = 1.5; /* This offsets the labels a little bit */
				foreach a of numlist 1(1)`nAge' {;
					foreach Spec in `SpecGraphList' {;
						local SpecXLabel `"`SpecXLabel' `s' "``Spec'Label'" "'; /* " */
						local s = `s' + 2;
					};
					local s = `s' + 1; 
				};
				 di `"`SpecXLabel'"'; /* " */
					
			/*--------------*/
			/* * * Plot * * */
			/*--------------*/
			
			twoway `PlotStatement'
					(rcap upper lower x, color(gs0) lcolor(black) lwidth(.1))
					(scatter sig2 x if ged == 1, msymbol(o) mlcolor(black) mfcolor(white) msize(huge) )
					(scatter sig1 x if ged == 1, msymbol(o) mcolor(black) msize(large) )
					(scatter sig1 x if hsg == 1, msymbol(d) mcolor(black) msize(large) )
					(scatter min x,   msymbol(o) mcolor(none) msize(tiny) )
					(scatter max x,   msymbol(o) mcolor(none) msize(tiny) ),
					legend( rows(2) order( 1 "GED" 4 "p<0.05 (GED vs.HSG)" 5 "p<0.05 (GED vs.Drop)" 2 "HSG" 6 "p<0.05 (HSG vs.Drop)" 3 "S.E."  ) symxsize(*.75)) yline(0, lcolor(black) ) xtitle("", size(small))
					xlabel( `AgeXLabel'  ,  angle(0)  labsize(medlarge) noticks custom labgap(*15) )
					xlabel( `SpecXLabel' ,  angle(-45)  labsize(medsmall) noticks custom add )
					title("") scheme(s1manual) xtitle("", size(small)) yline(0, lcolor(black)) ;				
			graph export "`FIG_DIR'/RawAbilBgNLSY79Selec_`DepVar'_`Pse'_`GenderVar'_`RaceVar'.png", replace;
			graph export "`FIG_DIR'/RawAbilBgNLSY79Selec_`DepVar'_`Pse'_`GenderVar'_`RaceVar'.eps", replace fontface(times) logo(off);
			
			#delimit cr;
		
	} /* End of Loop Across Dependent Variables */
	} /* End of Loop Across Genders */
	} /* End of Loop Across Race */
	} /* End of Loop Across Pse Runs */


} /* End Graph Maker */




