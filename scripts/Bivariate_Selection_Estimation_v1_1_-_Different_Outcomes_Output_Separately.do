*===============================================================================
*
*    Generate Graphs Comparing Estimation for NLSY79 Including Selection Correction for Labor Force Participation
*
*    By Nicholas Mader 12/8/2011
*
*===============================================================================


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

	local RunEstimation = 1
	local RunGraphing   = 0

	local INPUT_DIR  "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Intermediate Data/"
	local ROOT_DIR   "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Cross-Sectional Estimates"
	cd "`ROOT_DIR'"
	local FIG_DIR    "`ROOT_DIR'/output/figures"
	local TAB_DIR    "`ROOT_DIR'/output/tables"
	local MAT_DIR    "`ROOT_DIR'/output/matrices"

	local NLSY79data "`INPUT_DIR'/NLSY79 v8_5 - Partially Rebuilt Data to 2008 - Added Var Construction - Has Factors.dta"
	use "`NLSY79data'", clear
	
	foreach v of varlist asset* {
		replace `v' = `v'/10000
	}

/*============================*/
/*                            */
/* * * SET RUN PARAMETERS * * */
/*                            */
/*============================*/
	
	/* * * Run Parameters * * */
	
		local OutCVarList     = " LvlW " /*  Empl LvlHrW LvlH */
			local LvlWLabel   = "Annual Earnings"
			local EmplLabel   = "Probability of Employment"
			local LvlHrWLabel = "Hourly Wages"
			local LvlHLabel   = "Annual Hours Worked"
		local RaceVarList   = " all white black hisp " /*  */
		local GenderVarList = " female male " /*   */
		local AgeGroupList  = " A2024 A2529 A3034 A3539 " /*  */
			local A2024Label = "Age 20 to 24" 
			local A2529Label = "Age 25 to 29"
			local A3034Label = "Age 30 to 34"
			local A3539Label = "Age 35 to 39"
		local SpecEstList   = "AbilBg Cond Sel1 Sel2" /* Raw  */
		local SpecGraphList = "AbilBg Cond Sel1 Sel2" /* Raw  */
			local AbilbgLabel = "Abil+BG"
			local Cond = "Cond'l" 
			local Sel1 = "Univar Sel"
			local Sel2 = "Bivar Sel"
		local PseRunList = "AC" /* NC */

	/* * * Controls for the Dual Selection Procedure * * */
	
		global SwitchVar   ged_now
		local  SwitchCtrls afqt_pre_std mhgc famstatus_broken famstatus_hybrid south_age14 urban_age14 age ///
				married_now num_children_in_hh baby_in_hh toddler_in_hh SpouseWage unemp1-unemp4 lowage1-lowage4   preg_drop gov_affil_yet /*  */
			
		/* global EmplVar     LFP .... This assignment will be determined within OutCVar loop. Selection on LFP goes with LvlW and Empl, and Empl goes with LvlHrW and LvlH. */
		local  EmplCtrls0  afqt_pre_std mhgc famstatus_broken famstatus_hybrid south_age14 urban_age14 age ///
				married_now num_children_in_hh baby_in_hh toddler_in_hh SpouseWage unemp1-unemp4 lowage1-lowage4  asset_tot*_ti /*  asset_prop*_ti asset_fin*_ti */
		local  EmplCtrls1  `EmplCtrls0'
		
		/* local  OutCVar     LvlW ..... This assignment will be determined within the OutCVar loop. */
		local  OutCCtrls   afqt_pre_std mhgc famstatus_broken famstatus_hybrid south_age14 urban_age14 age dregion2 dregion3 dregion4 _Iyear* black hisp
		local  RawCtrls    age dregion2 dregion3 dregion4 _Iyear* black hisp
		
		local SampleRestr     = "jail_yet!=1 & enrolled!=1 "  /* & (LvlH <=4000 | LvlH ==.) & (LvlW<300000 | LvlW==.) & LvlHrW>=3 & LvlHrW<=200 */
		local PseAC_Restr     = ""
		local PseNC_Restr     = "& (touchcoll_ever == 0)"
		local PseOC_Restr     = "& (touchcoll_ever == 1)"

	/* * * Get Dual Selection Estimation Code Loaded into Memory * * */
	
		run "./scripts/Dual Selection Method Code - 1 - Run Dual Regime Selection Method - v2 - Macros converted to globals.do"
		run "./scripts/Dual Selection Method Code - 2 - Run Quadrature - v2 - Macros converted to globals.do"

		/* Set up program to save results from the outcome equation */
		cap program drop RunTests
		program RunTests
			syntax anything [, row(string) spec(string) age(string) depvar(string) matname(string)]
			matrix `matname'[`row', 1] = _b[ged_now]
			matrix `matname'[`row', 2] = _se[ged_now]
			test ged_now = 0
			if r(p) < 0.05 matrix `matname'[`row', 3] = _b[ged_now]
			matrix `matname'[`row', colnumb(`matname', "`spec'")]   = 1
			matrix `matname'[`row', colnumb(`matname', "`age'")]    = 1
			matrix `matname'[`row', colnumb(`matname', "`depvar'")] = 1
		end
	
	/*--------------------------------*/
	/* RUN LOOPS ACROSSS DEMOGRAPHICS */
	/*--------------------------------*/
	
	
	foreach GenderVar in `GenderVarList' {
		local Gx = "G" + substr("`GenderVar'",1,1)
	foreach RaceVar in `RaceVarList' {
		local Rx = "R" + substr("`RaceVar'",1,1)
		if "`RaceVar'" == "all" local SubSampleRestr = ""
		if "`RaceVar'" != "all" local SubSampleRestr = " & CrossSectSample == 1 "
	foreach Pse in `PseRunList' {
	
		di "Now Running GenderVar = `GenderVar', RaceVar = `RaceVar', Pse = `Pse'."

	#delimit;
	foreach OutCVar in `OutCVarList' {;
	
		if "`OutCVar'" == "LvlW" | "`OutCVar'" == "Empl"   global EmplVar LFP;
		if "`OutCVar'" == "LvlH" | "`OutCVar'" == "LvlHrW" global EmplVar Empl;
		
		/*------------------------------------*/
		/* Generate Matrix for Saving Results */
		/*------------------------------------*/
		
		local nAges    = wordcount("`AgeGroupList'");
		local nSpecs   = wordcount("`SpecEstList'");
		local nSlotsPerGraph = (`nAges'+1)*(`nSpecs')-1; /* The "+1" allows room for a blank space after each age group. "nSpecs" indicates the size of each group of plotted estimates. */
		local RegColNames beta se sig `SpecEstList' `AgeGroupList';
		local RegMat SelectEsts_`OutCVar'_`Gx'`Rx'`Pse';
		matrix `RegMat' = J(`nSlotsPerGraph', wordcount("`RegColNames'"), .);
		matrix colnames `RegMat' = `RegColNames';
		
		local r = 1;
		eststo clear;
	
		foreach Axxyy in `AgeGroupList' {;
		
			preserve;
			keep if `SampleRestr' `SubSampleRestr' & `GenderVar' == 1 & `Axxyy' == 1 & (dropout_now == 1 | ged_now == 1) `Pse`Pse'_Restr';
			keep id year $SwitchVar `SwitchCtrls' `EmplVar' `EmplCtrls0' `EmplCtrls1' `OutCVar' `OutCCtrls';

			/*-----------------------------*/
			/***   BIVARIATE SELECTION   ***/
			/*-----------------------------*/
			
			if strpos("`SpecEstList'", "Sel2") > 0 {;
				
				di "Contruction Control Function Values for Bivariate Selection Correction";
					
				/*--------------------------*/
				/* Obtain Initial Estimates */
				/*--------------------------*/
					
				local start_time = c(current_time);
				qui count;
				local nEst = r(N);
					
					probit $SwitchVar `SwitchCtrls', robust;
						cap drop zg;
						predict zg, xb;
					probit $EmplVar `EmplCtrls0' if $SwitchVar == 0, robust;
						matrix b2_0 = e(b);
					qui probit $EmplVar `EmplCtrls1' if $SwitchVar == 1, robust;
						matrix b2_1 = e(b);
					corr $EmplVar $SwitchVar;
						local logitrho = logit((r(rho)+1)/2);
					matrix b2_init = b2_0, b2_1, scalar(`logitrho'), scalar(`logitrho');
				
				/*-----------------------------------*/
				/* Run Dual Hiring Regime Estimation */
				/*-----------------------------------*/
				
					global PredZG zg;
					global Ind_y1 $SwitchVar;
					ml model lf CondProb_2stage (h0: $EmplVar = `EmplCtrls0') (h1: $EmplVar = `EmplCtrls1') /lrho0 /lrho1, robust
						diparm(lrho0, function(invlogit(@)*2-1) derivative(( invlogit(@)*(1-invlogit(@)) )*2) label("Non GED Rho"))
						diparm(lrho1, function(invlogit(@)*2-1) derivative(( invlogit(@)*(1-invlogit(@)) )*2) label("Has GED Rho"));
					ml init b2_init, copy;
					ml maximize;
					
				di "456";
					
					cap drop zf0 zf1;
					predict zf0, equation(h0) xb;
					predict zf1, equation(h1) xb;
				
					nlcom invlogit(_b[lrho0:_cons])*2-1;
						global rho0    = el(r(b),1,1);
						global rho0_se = sqrt(el(r(V),1,1));
					nlcom invlogit(_b[lrho1:_cons])*2-1;
						global rho1    = el(r(b),1,1);
						global rho1_se = sqrt(el(r(V),1,1));


			/*------------------------------------------*/
			/***  ESTIMATE CONTROL FUNCTIONS VALUES   ***/
			/*------------------------------------------*/

				/* * * Obtain values for the control functions associated with the dual selection method * * */

					RunDualSelectionQuad;
					
			local end_time = c(current_time);
			di "Bivariate selection correction method for `OutCVar' and `Axxyy' started at `start_time' and ended at `end_time', having run `nEst' cases.";		
			
		}; /* End of Construction of Bivariate Control Variables */
		
			/* Run Analyses if Called For */
			
			if strpos("`SpecEstList'","Raw") > 0 {;
				eststo Raw_`OutCVar'_`Axxyy' : reg     `OutCVar' $SwitchVar `RawCtrls',                   vce(cluster id);
				RunTests hello, row(`r') spec("Raw") age("`Axxyy'") depvar("`OutCVar'") matname("`RegMat'"); local r = `r' + 1;
			};
				
			if strpos("`SpecEstList'","AbilBg") > 0 {;
				eststo Reg_`OutCVar'_`Axxyy' : reg     `OutCVar' $SwitchVar `OutCCtrls',                  vce(cluster id);
				RunTests hello, row(`r') spec("AbilBg") age("`Axxyy'") depvar("`OutCVar'") matname("`RegMat'"); local r = `r' + 1;
			};
				
			if strpos("`SpecEstList'","Cond") > 0 {;
				eststo Cond_`OutCVar'_`Axxyy': reg     `OutCVar' $SwitchVar `OutCCtrls' if $EmplVar == 1, vce(cluster id);
				RunTests hello, row(`r') spec("Cond") age("`Axxyy'") depvar("`OutCVar'") matname("`RegMat'"); local r = `r' + 1;
			};
			if strpos("`SpecEstList'","Sel1") > 0 {;
				eststo Sel1_`OutCVar'_`Axxyy': heckman `OutCVar' $SwitchVar `OutCCtrls', select($EmplVar = `EmplCtrls1') twostep;
				RunTests hello, row(`r') spec("Sel1") age("`Axxyy'") depvar("`OutCVar'") matname("`RegMat'"); local r = `r' + 1;
			};
			if strpos("`SpecEstList'","Sel2") > 0 {;
				eststo Sel2_`OutCVar'_`Axxyy': reg     `OutCVar' $SwitchVar `OutCCtrls' EV_Dual_G0 EV_Dual_G1 EM0_Dual_G0 EM1_Dual_G1, vce(cluster id);
				RunTests hello, row(`r') spec("Sel2") age("`Axxyy'") depvar("`OutCVar'") matname("`RegMat'"); local r = `r' + 1;
			};
				
			local r = `r' + 1;
				
			* pause;
			
		}; /* End of Loop Across Age */
		
		restore
	
	}; /* End of Loop Across Outcome Variables */
	
		#delimit cr;
	
		matsave `RegMat',  path("`MAT_DIR'") sav replace
		esttab using "`TAB_DIR'/NLSY79 - Bivariate Selection and Comparisons - `GenderVar', `RaceVar', `Pse'.csv", compress nostar replace csv nonumbers mtitles se
	
	} /* End of Loop Across Pse */
	} /* End of Loop Across Race */
	} /* End of Loop Across Gender */

	
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
	foreach OutCVar in `OutCVarList' {
		
		/* Identify the Minimum and Maximum Estimates Across Age and Specification (but within Gender, given how disparate the estimates are)
		 
		local max_y = 0
		local min_y = 0
		
		foreach AgeGrp in `AgeGroupList' {
			foreach Spec in `SpecGraphList' {
				foreach edu in ged hsg {
					local reg_mat SelectEsts_`Gx'`Rx'`Pse'
					matload `reg_mat', path("`MAT_DIR'") overwrite sav
					keep if `OutCVar' == 1
					local new_beta = `reg_mat'[rownumb(`reg_mat',"`data'`edu'"), colnumb(`reg_mat',"beta")]
					local new_se   = `reg_mat'[rownumb(`reg_mat',"`data'`edu'"), colnumb(`reg_mat',"se")]
					if `new_beta'+`new_se'>`max_y' local max_y = `new_beta'+`new_se'
					if `new_beta'-`new_se'<`min_y' local min_y = `new_beta'-`new_se'
				} /* End of Loop Across Education */
			} /* End of Loop Across specifications */
		} /* End of Loop Across Ages */
		
		#delimit;
		matload `reg_mat', path("`MAT_DIR'") overwrite sav missing(-999);
		keep if "`OutCVar'" == 1;
		foreach v of * {; replace `v' = . if `v' == -999; };

		/* Generate New Variables for Plotting */
			gen x       = _n;
			/* gen max_val = `max_val';
			gen min_val = `min_val'; */
			local nAge  = wordcount("`AgeGroupList'");
			local nSpec = wordcount("`SpecGraphList'");
		
		/* Series Plot Statement */
			local PlotStatement;
			if nSpec == 3 {; local fcolor1 white; local fcolor2 gray; local fcolor3 black; };
			if nSpec == 4 {; local fcolor1 white; local fcolor2 gs4;  local fcolor3 gs12;  local fcolor4 black; };
			local i = 1;
			foreach Spec in `SpecGraphList' {;
				local PlotStatement `PlotStatement' (bar beta x if `Spec' == 1, lcolor(black) fcolor(`fcolor`i'') lwidth(0.25));
				local i = `i' + 1;
			};
		
		/* Legend Order Statement */
			local LegendOrder;
			local s = 1;
			foreach Spec in `SpecGraphList' {;
				local LegendOrder `"`LegendOrder' `s' "``Spec'Label'" "'; /* " */
				local s = `s' + 1;
			};
			local sp1 = `s' + 1;
			local LegendOrder `" `LegendOrder' `s' "+/- 1 S.E." `sp1' "5% Sig (GED vs. DO)" "';  /* " */
			* di `"`LegendOrder'"';
			
		/* x Label for Ages */
			local xLabel;
			local a = `nSpec'/2; /* This offsets the labels to lay centered below the data series */
			foreach AgeGroup in `AgeGroupList' {;
				local xLabel `" `xLabel' `a' "``AgeGroup'Label'" "'; /* " */
				local a = `a' + `nSpec' + 1;
			};
			di `"`xLabel'"'; /* " */

		
		twoway `PlotStatement'
			(rcap upper lower x , color(gs0) lcolor(black) lwidth(.1))
			(sc sig, msymbol(o) mcolor(black) msize(large) )
			/* (sc min model,   msymbol(o) mcolor(none) msize(tiny) )
			(sc max model,   msymbol(o) mcolor(none) msize(tiny) ) */
		, legend( rows(2) order( `LegendOrder') ) xlabel( `xLabel' ,  angle(0) noticks custom)
		title("")  scheme(s1manual)  yline(0, lcolor(black)) xtitle("", size(small)) ``OutC'Label'; 
		
		graph export "`FIG_DIR'/RawAbilBg`Pse'NLSY79Selec_`DepVar'_`GenderVar'_`RaceVar'.png", replace;
		graph export "`FIG_DIR'/RawAbilBg`Pse'NLSY79Selec_`DepVar'_`GenderVar'_`RaceVar'.eps", replace logo(off);
		
		#delimit cr;
	} /* End of Loop Across Genders */
	} /* End of Loop Across Race */
	} /* End of Loop Across Dependent Variables */
	} /* End of Loop Across Pse Runs */


} /* End Graph Maker */




