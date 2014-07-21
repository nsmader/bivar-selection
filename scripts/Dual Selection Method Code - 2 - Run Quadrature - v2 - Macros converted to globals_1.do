program RunDualSelectionQuad

	#delimit;
	local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
	local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 

	local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
	local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;

	cap drop EV_Dual EM_Dual EV_Marg EM_Marg;
	cap gen EV_Dual = .;  label var EV_Dual "Expected value of GED error, conditional on GED status and work.";
	cap gen EM_Dual = .;  label var EM_Dual "Expected value of work error, conditional on GED status and work.";
	
	cap gen EV_Marg = .;  label var EV_Marg "Expected value of GED error, conditional on GED status only.";
	cap gen EM_Marg = .;  label var EM_Marg "Expected value of work error, conditional on employment status only.";
	
	qui count;
	local N = r(N);
	
	local start_time = c(current_time);
	
	forvalues i = 1(1)`N' {;
	
		if mod(`i',floor(`N'/50)) == 0 di "Integration is " `i'/`N'*100 "% complete."; 

		foreach ExpectVar in V M {; /* This integral procedure will calculate the expected M# for the appropriate A=0 or A=1 outcome */
			
			local var = lower("`ExpectVar'");
		
			/*----------------------------*/
			/* Set Bounds for Integration */
			/*----------------------------*/
			
				local SwitchVal = $SwitchVar[`i'];
				local EmplVal   = ${EmplVar}[`i'];
				local rho_i     = ${rho`SwitchVal'};
									
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
	
	di "Integration took from `start_time' to " c(current_time);
	
	cap drop EM0_Dual EM1_Dual EM0_Marg EM1_Marg;
	gen EM0_Dual = 0; gen EM1_Dual = 0; gen EM0_Marg = 0; gen EM1_Marg = 0;
	
	label var EV_Dual "Expected value of GED error, conditional on GED and work outcome.";
	label var EV_Marg "Expected value of GED error, conditional on GED outcome.";
	
	replace EM0_Dual = EM_Dual if $SwitchVar == 0 & $EmplVar == 1; label var EM1_Dual "Expected value of work error, conditional on GED and Work.";
	replace EM1_Dual = EM_Dual if $SwitchVar == 1 & $EmplVar == 1; label var EM0_Dual "Expected value of work error, conditional on no GED and Work.";
	
	replace EM1_Marg = EM_Marg if $EmplVar == 1; label var EM0_Marg "Expected value of work error, conditional on Working.";
	replace EM0_Marg = EM_Marg if $EmplVar == 0; label var EM0_Marg "Expected value of work error, conditional on Not Working.";
	
	/* Make sure that all of the control function variables only take on values for the case that they apply to. */
	cap drop E*_Dual_G*; 
	gen EV_Dual_G0  = EV_Dual * ($SwitchVar == 0) * ($EmplVar == 1);
	gen EV_Dual_G1  = EV_Dual * ($SwitchVar == 1) * ($EmplVar == 1);
	gen EM0_Dual_G0 = EM_Dual * ($SwitchVar == 0) * ($EmplVar == 1);
	gen EM1_Dual_G1 = EM_Dual * ($SwitchVar == 1) * ($EmplVar == 1);
	#delimit cr;
	
end
