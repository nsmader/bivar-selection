/*-------------------------------------------------------*/
/*** ESTIMATION FOR SELECTION INTO GED AND EMPLOYMENT  ***/
/*-------------------------------------------------------*/

	clear all
	set mem 1g
	set more off
	pause on
		
	local SwitchVar   ged_now
	local SwitchCtrls afqt_pre_std famstatus_broken famstatus_hybrid south_age14 urban_age14 married_now age age2 baby_in_hh toddler_in_hh SpouseWage unemp1-unemp4 preg_drop govt_train
	local EmplVar     Empl
	local EmplCtrls0  afqt_pre_std famstatus_broken famstatus_hybrid south_age14 urban_age14 married_now age age2 baby_in_hh toddler_in_hh SpouseWage assets_fin assets_prop unemp1-unemp4
	local EmplCtrls1  `EmplCtrls0'
	local OutCVar     LvlHrW 
	local OutCCtrls   afqt_pre_std famstatus_broken famstatus_hybrid south_age14 urban_age14 ActExp ActExp2
		
		
	/* MODEL: 
		A = zg - V,             (this is the "switching" outcome; in our application, A represents GED certification)
		B = zf0 - M0 if A = 0
		B = zf1 - M1 if A = 1   (this is the sample selection outcome; in our application, B represents the employment outcome)
		Y = xb + d*A + U        (this is the observed outcome of interest, observed only if B=1)
		U, V, M0, and M1 are distributed multivariate normal.
	*/
	
	/* TO DO:
		o Is there an easier way to directly extract the transformed rho0 and rho1 values that are currently output by _diparm?
		o Use own calculation of bivariate normal cumulative densities instead of monkeying with predefined function
	*/

/*-------------------------------*/
/***   DUAL REGIME SELECTION   ***/
/*-------------------------------*/

	cap program drop CondProb_2stage
	program CondProb_2stage
	
		args lf zf0 zf1 logitrho0 logitrho1
			/* Notes:
				(1) logitrho terms equal logit((rho+1)/2) so that we represent the parameter on the real line;
				(2) zf0 and zf1 are observed predictor terms for probability of B when, respectively, A equals 0 and A equals 1;
			*/
			
		/* Note: In this function, choice 1 represents the decision to GED certify. This is not modeled here, but rather is controlled by a previous stage's estimation,
			and the predicted index value zg is stored in a data variable. Choice 2 is the choice to work or not, conditional on choice 1. The parameters governing choice 2
			conditional on GED=0 and GED=1 are predicted here, and are represented in the indices "zf0" and "zf1" (short for "z phi 0" and "z phi 1"). */  
		
		tempvar rho0 rho1 p00 p01 p10 p11
		
		qui gen double `rho0' = invlogit(`logitrho0')*2 - 1
		qui gen double `rho1' = invlogit(`logitrho1')*2 - 1
		
		
		/* Because the binormal(.,.) function only calculates cumulative density below given threshold values, we calculate upper-tail densities using the following
			properties of the normal distribution (due to symmetry):
				Pr(V<v, M>m, rho(V,M)) = Pr(V<v, -M<-m, rho(V,-M)) = binormal(v,-m,-rho(V,M)) 
					and
				Pr(V>v, M>m, rho(V,M)) = Pr(-V<-v, -M<-m, rho(-V,-M)) = binormal(-v,-m,rho(V,M))  */

		qui gen double `p00'   = binormal(-$PredZG, -`zf0',  `rho0')
		qui gen double `p01'   = binormal(-$PredZG,  `zf0', -`rho0')
		qui gen double `p10'   = binormal( $PredZG, -`zf1', -`rho1')
		qui gen double `p11'   = binormal( $PredZG,  `zf1',  `rho1')
		
		qui replace `lf' = (1-$Ind_y1) * (1-$ML_y1) * ln(`p00') ///
						 + (1-$Ind_y1) * (  $ML_y1) * ln(`p01') ///
						 + (  $Ind_y1) * (1-$ML_y2) * ln(`p10') ///
						 + (  $Ind_y1) * (  $ML_y2) * ln(`p11')
	end
		

	/*--------------------------*/
	/* Obtain Initial Estimates */
	/*--------------------------*/
		
		qui probit `SwitchVar' `SwitchCtrls', robust
		cap drop zg
		predict zg, xb
		
		qui probit `EmplVar' `EmplCtrls0' if `SwitchVar' == 0, robust
		matrix b2_0 = e(b)
		
		qui probit `EmplVar' `EmplCtrls1' if `SwitchVar' == 1, robust
		matrix b2_1 = e(b)
		
		corr `EmplVar' `SwitchVar'
		local logitrho = logit((r(rho)+1)/2)
		
		matrix b2_init = b2_0, b2_1, scalar(`logitrho'), scalar(`logitrho')
	
	
	/*-----------------------------------*/
	/* Run Dual Hiring Regime Estimation */
	/*-----------------------------------*/
	
		global PredZG zg
		global Ind_y1 `SwitchVar'
		ml model lf CondProb_2stage (h0: `EmplVar' = `EmplCtrls0') (h1: `EmplVar' = `EmplCtrls1') /lrho0 /lrho1, robust ///
			diparm(lrho0, function(invlogit(@)*2-1) derivative(( invlogit(@)*(1-invlogit(@)) )*2) label("Non GED Rho")) ///
			diparm(lrho1, function(invlogit(@)*2-1) derivative(( invlogit(@)*(1-invlogit(@)) )*2) label("Has GED Rho"))
		ml init b2_init, copy
		ml maximize
		
		/* Save Results */
		
			cap drop zf0 zf1
			predict zf0, equation(h0) xb
			predict zf1, equation(h1) xb 
		
			nlcom invlogit(_b[lrho0:_cons])*2-1
				local rho0    = el(r(b),1,1)
				local rho0_se = sqrt(el(r(V),1,1))
				
			nlcom invlogit(_b[lrho1:_cons])*2-1
				local rho1    = el(r(b),1,1)
				local rho1_se = sqrt(el(r(V),1,1))			

				
/*------------------------------------------*/
/*------------------------------------------*/
/***  ESTIMATE CONTROL FUNCTIONS VALUES   ***/
/*------------------------------------------*/
/*------------------------------------------*/

	#delimit;
	local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
	local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 

	local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
	local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;
		
	cap gen EV_Dual = .;  label var EV_Dual "Expected value of GED error, conditional on GED status and work.";
	cap gen EM_Dual = .;  label var EM_Dual "Expected value of work error, conditional on GED status and work.";
	
	cap gen EV_Marg = .;  label var EV_Marg "Expected value of GED error, conditional on GED status only.";
	cap gen EM_Marg = .;  label var EM_Marg "Expected value of work error, conditional on employment status only.";
	
	preserve;
	keep if `SampleRestr';
	qui count;
	local N = r(N);
	
	* local start_time = c(current_time);
	
	forvalues i = 1(1)`N' {;
	
		* if mod(`i',floor(`N'/50)) == 0 di "Integration is " `i'/`N'*100 "% complete."; 

		foreach ExpectVar in V M {; /* This integral procedure will calculate the expected M# for the appropriate A=0 or A=1 outcome */
			
			local var = lower("`ExpectVar'");
		
			/*----------------------------*/
			/* Set Bounds for Integration */
			/*----------------------------*/
			
				local SwitchVal = `SwitchVar'[`i'];
				local EmplVal   = `EmplVar'[`i'];
				local rho_i     = `rho`SwitchVal'';
									
				local V_TruncVal = zg[`i'];            local Upper_V_TruncInd = `SwitchVal'; /* Note: Upper_V_TruncInd = GED corresponds to GED=1 <=> zg-v>0 <=> v<zg which is indeed an indication of upper truncation at zg. */ 
				local M_TruncVal = zf`SwitchVal'[`i']; local Upper_M_TruncInd = `EmplVal';
				
				/* Set integration bounds for our cdf transformations -- u_v and u_m -- of their corresponding normal variables */ 
					/* Note that the "@_mult" terms above will be used in our trick in calling the binormal(.,.) function that is described above */
				
				if `Upper_V_TruncInd' == 1 {; local a_v = 0.0;                  local b_v = normal(`V_TruncVal'); local v_mult =  1; };
				if `Upper_V_TruncInd' == 0 {; local a_v = normal(`V_TruncVal'); local b_v = 1.0;                  local v_mult = -1; };
				if `Upper_M_TruncInd' == 1 {; local a_m = 0.0;                  local b_m = normal(`M_TruncVal'); local m_mult =  1; };
				if `Upper_M_TruncInd' == 0 {; local a_m = normal(`M_TruncVal'); local b_m = 1.0;                  local m_mult = -1; };
					
				/* Set up Variable Transformation for Gaussian Sum. */
			
				local t1_v = ((`b_v')-(`a_v'))/2; local t2_v = ((`b_v')+(`a_v'))/2;
				local t1_m = ((`b_m')-(`a_m'))/2; local t2_m = ((`b_m')+(`a_m'))/2;

			/*---------------------*/			
			/* Run the Integration */
			/*---------------------*/

				local QSum = 0;
				local QSum_Marg = 0;
				
				forvalues p1 = 1(1)10 {;
				
					/* Calculate marginal expectations of each variable, i.e. E[V|V>a] and E[M|M>b] */
					local ExpectVar_i_marg = invnormal((`t1_`var'')*(`x`p1'') + (`t2_`var''));
					local QSum_Marg = `QSum_Marg' + `w`p1''*( (0.5)*`ExpectVar_i_marg' );
				
					forvalues p2 = 1(1)10 {;	
				
						local V_i = invnormal((`t1_v')*(`x`p1'') + (`t2_v'));
						local M_i = invnormal((`t1_m')*(`x`p2'') + (`t2_m')); /* Note that this "M" variable generically stands for either "M0" or "M1", based on what the individual's GED status is. */
						
						local binormpdf_num = exp( (-1) * (((`V_i')^2) + ((`M_i')^2) - (2*(`rho_i')*(`V_i')*(`M_i'))) / (2*(1-(`rho_i')^2)) );  /* Numerator of the bivariate normal */
						local binormpdf_den = 2*_pi*sqrt(1-(`rho_i')^2);                                                                            /* Denominator of the bivariate normal */
						local binormpdf = `binormpdf_num'/`binormpdf_den';
						
						local Summand =  (``ExpectVar'_i'*`binormpdf') / (normalden(`V_i')*normalden(`M_i'));
						local QSum = `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' );
						
					};
				};
				
				local binorm_prob = binormal((`v_mult')*(`V_TruncVal'), (`m_mult')*(`M_TruncVal'), (`v_mult')*(`m_mult')*(`rho_i'));
				local QSum = (`t1_v')*(`t1_m')*(`QSum')/`binorm_prob';
				
				qui replace E`ExpectVar'_Dual = `QSum'      if _n == `i';
				qui replace E`ExpectVar'_Marg = `QSum_Marg' if _n == `i';
			
		}; /* End of Loop Calculating Expected Value of Both V and M */
		
	}; /* End of Loop Across Individuals */
	
	* di "Integration took from `start_time' to " c(current_time);
	
	cap drop EM0_Dual EM1_Dual EM0_Marg EM1_Marg;
	gen EM0_Dual = 0; gen EM1_Dual = 0; gen EM0_Marg = 0; gen EM1_Marg = 0;
	
	label var EV_Dual "Expected value of GED error, conditional on GED and work outcome.";
	label var EV_Marg "Expected value of GED error, conditional on GED outcome.";
	
	replace EM0_Dual = EM_Dual if `SwitchVar' == 0 & `EmplVar' == 1; label var EM1_Dual "Expected value of work error, conditional on GED and Work.";
	replace EM1_Dual = EM_Dual if `SwitchVar' == 1 & `EmplVar' == 1; label var EM0_Dual "Expected value of work error, conditional on no GED and Work.";
	
	replace EM1_Marg = EM_Marg if `EmplVar' == 1; label var EM0_Marg "Expected value of work error, conditional on Working.";
	replace EM0_Marg = EM_Marg if `EmplVar' == 0; label var EM0_Marg "Expected value of work error, conditional on Not Working.";

/*----------------------------------------------------------------*/		
/* Run Outcome Estimation With Selection and Endogenous Regressor */
/*----------------------------------------------------------------*/			

	/* Make sure that all of the control function variables only take on values for the case that they apply to. */
	cap drop E*_Dual_G*; 
	gen EV_Dual_G0  = EV_Dual * (`SwitchVar' == 0) * (`EmplVar' == 1);
	gen EV_Dual_G1  = EV_Dual * (`SwitchVar' == 1) * (`EmplVar' == 1);
	gen EM0_Dual_G0 = EM_Dual * (`SwitchVar' == 0) * (`EmplVar' == 1);
	gen EM1_Dual_G1 = EM_Dual * (`SwitchVar' == 1) * (`EmplVar' == 1);
	
	reg `OutCVar' `OutCCtrls' `SwitchVar' EV_Dual_G0 EV_Dual_G1  EM0_Dual_G0  EM1_Dual_G1, vce(cluster id);
	
	restore;
			
#delimit cr;
