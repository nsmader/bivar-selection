/*--------------------------------------------*/
/* TEST CODE FOR CONTROL FUNCTION CORRECTIONS */
/*--------------------------------------------*/

/*-----------------------*/
/* GENERATE TOY DATA SET */
/*-----------------------*/

	clear all
	set mem 1g
	set more off
	pause on

	/* Generate Variables */
	
		set obs 10000
		set seed 1092011
	
		gen x1 = runiform()
		gen x2 = runiform()
		gen x3 = runiform()
		gen x4 = runiform()
	
	/* Generate Errors */
	
		/* Can find notes for this at http://www.stata.com/support/faqs/stat/mvnorm.htmlhttp://www.stata.com/support/faqs/stat/mvnorm.html */
	
		gen c1 = invnorm(runiform())
		gen c2 = invnorm(runiform())
		gen c3 = invnorm(runiform())
		gen c4 = invnorm(runiform())
		
		matrix Sig_u = (1.0, 0.5, 0.5, 0.3 \ ///
		                0.5, 1.0, 0.8, 0.7 \ ///
		                0.5, 0.8, 1.0, 0.7 \ ///		                
		                0.3, 0.7, 0.7, 1.0 )
		             
		matrix Chol_Sig_u = cholesky(Sig_u)
		
		local rows = rowsof(Chol_Sig_u)
		
		foreach r of numlist 1(1)`rows' {
			matrix CSu`r' = Chol_Sig_u[`r',1...]	
			matrix score u`r' = CSu`r' /* This works because the default matrix names are c# and that's what we've used as our variable names for orthogonal errors */
		} 
 
	/* Generate Outcomes */

		gen     y1  = (-1.5 + 3.0*x1 + u1 > 0)             /* This is whether an individual gets a GED */
		
		gen     y2_1regime  = ( 0.3 - 2.0*x2 + u2 > 0)
		gen     y2_2regime  = ( 0.3 - 2.0*x2 + u2 > 0) if y1 == 0  /* This is whether an individual works or not, conditional on not having a GED */
		replace y2_2regime  = ( 0.4 - 1.0*x3 + u3 > 0) if y1 == 1  /* This is whether an individual works or not, conditional on     having a GED */ 
		
		gen     y3_1regime = .
		gen     y3_2regime = .
		replace y3_1regime = 0.0 + 3.0*x4 + 5.0*y1 + u4 if y2_1regime == 1
		replace y3_2regime = 0.0 + 3.0*x4 + 5.0*y1 + u4 if y2_2regime == 1
		
		/* Generate an outcome that only features endogeneity with respect to based on y1 */ 
		gen     y4_endogtreat = .
		replace y4_endogtreat = 0.0 + 3.0*x4 + 5.0*y1 + u4
		
		/* Generate an outcome that only features sample selection with respect to y2 */
		
		gen     y4_sel_1regime = .
		replace y4_sel_1regime = 0.0 + 3.0*x4 + u4 if y2_1regime == 1
		

/*---------------------------------------------------*/
/*---------------------------------------------------*/
/***   RUN DOUBLE SELECTION CORRECION ESTIMATION   ***/
/*---------------------------------------------------*/
/*---------------------------------------------------*/

	/* See the files "Sandbox Code for Running MLE.do" and "Sandbox Code for Numerical Integration.do" for examples of simpler estimation methods using STATA's MLE module, and
		for examples of numerical integration using both Monte Carlo and Gaussian Quadrature methods. */

	/*-------------------*/
	/* Define Estimation */
	/*-------------------*/
	
		cap program drop CondProb_2stage
		program CondProb_2stage
		
			args lf zf0 zf1 logitrho0 logitrho1  /* Note: the logitrho terms equal logit((rho+1)/2) so that we represent the parameter on the real line */
				
			/* Note: In this function, choice 1 represents the decision to GED certify. This is not modeled here, but rather is controlled by a previous stage's estimation,
				and the predicted index value zg is stored in a data variable. Choice 2 is the choice to work or not, conditional on choice 1. The parameters governing choice 2
				conditional on GED=0 and GED=1 are predicted here, and are represented in the indices "zf0" and "zf1" (short for "z phi 0" and "z phi 1"). */  
			
			tempvar rho0 rho1 p00 p01 p10 p11
			
			qui gen double `rho0' = invlogit(`logitrho0')*2 - 1
			qui gen double `rho1' = invlogit(`logitrho1')*2 - 1
			
			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			/* Note: it is necessary to make creative calls of the binormal() cdf, because it only cumulates to an upper bound (and does not permit a call to calculate density of the upper tail).  
				I am currently making use of the symmetry of the distribution, where the area above our desired upper bound is equivalent to the area below the negative of that upper bound.
				This is effectively equivalent to examining properties of the negative of that random variable. When using that transformation, it is necessary to reverse the sign of the rho,
				since the sign of the covariance between the two random variables is now reversed. Put more simply: 
				
					Pr(A<a, B>b, rho(A,B)) = Pr(A<a, -B<-b, rho(A,-B)) = binormal(a,-b,-rho(A,B)) 
					
						and
						
					Pr(A>a, B>b, rho(A,B)) = Pr(-A<-a, -B<-b, rho(-A,-B)) = binormal(-a,-b,rho(A,B))	
				
				In our model, 
					
					p11 = Pr(GED=1, Empl=1|Z) = Pr(zg - V > 0, zf - M > 0) = Pr(V < zg, M < zf) = Phi_2(zg, zf(GED), rho(V,M)).  (This is motivated by Professor Heckman's notes on the two-dimensional selection
					
				The other probabilities are suitably modified according to the logic in the above comment. */
			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			
			
			qui gen double `p00'   = binormal(-$PredZG, -`zf0',  `rho0')
			qui gen double `p01'   = binormal(-$PredZG,  `zf0', -`rho0')
			qui gen double `p10'   = binormal( $PredZG, -`zf1', -`rho1')
			qui gen double `p11'   = binormal( $PredZG,  `zf1',  `rho1')

			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
			/* Note that "$Ind_y1" is an indicator value of whether someone has y1=1. This is set by the declaration of a global variable to reference a variable in the data. */
			/* The "$ML_y1" and "$ML_y2" terms are automatically generated by the MLE function to represent values of the variables declared as the dependent variables in the 1st and 2nd equations  
			   fed to the MLE. In this case, these are the values of y2 (under regime GED = 0) and y2 (under regime GED = 1). The fact that they are labeled as "...y1" and "...y2" by default is 
			   unfortunate coincidence. */
			/*--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/  
			
			qui replace `lf' = (1-$Ind_y1) * (1-$ML_y1) * ln(`p00') ///
			                 + (1-$Ind_y1) * (  $ML_y1) * ln(`p01') ///
			                 + (  $Ind_y1) * (1-$ML_y2) * ln(`p10') ///
			                 + (  $Ind_y1) * (  $ML_y2) * ln(`p11')
		end
		

	/*----------------------*/
	/* Initialize Estimates */
	/*----------------------*/
		
		qui probit y1 x1, robust
		cap drop zg
		predict zg, xb
		
		qui probit y2_2regime x2 if y1 == 0, robust
		matrix b2_0 = e(b)
		
		qui probit y2_2regime x3 if y1 == 1, robust
		matrix b2_1 = e(b)
		
		corr y1 y2_2regime
		local logitrho = logit((r(rho)+1)/2)
		
		matrix b2_init = b2_0, b2_1, scalar(`logitrho'), scalar(`logitrho')
	
	
	/*-----------------------------------*/
	/* Run Dual Hiring Regime Estimation */
	/*-----------------------------------*/
	
		global PredZG zg
		global Ind_y1 y1
		ml model lf CondProb_2stage (h0: y2_2regime = x2) (h1:y2_2regime = x3) /lrho0 /lrho1, robust ///
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
				
			di "rho0 = " `rho0' ", rho1 = " `rho1'
			
			cap drop rho0 rho1
			cap gen rho0 = `rho0'
			cap gen rho1 = `rho1'
			

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
			
		cap gen EV_CtrlFn = .;  label var EV_CtrlFn  "Expected value of GED error, conditional on GED status and work.";
		cap gen EM_CtrlFn = .;  label var EM_CtrlFn "Expected value of work error, conditional on GED status and work.";
		
		cap gen EV_MargExpect = .;  label var EV_MargExpect "Expected value of GED error, conditional on GED status only.";
		cap gen EM_MargExpect = .;  label var EM_MargExpect "Expected value of work error, conditional on employment status only.";
		
		qui count;
		local N = r(N);
		
		local start_time = c(current_time);
		
		forvalues i = 1(1)`N' {;
		
			if mod(`i',floor(`N'/50)) == 0 di "Integration is " `i'/`N'*100 "% complete."; 

			foreach ExpectVar in V M {;
				
				local var = lower("`ExpectVar'");
			
				/* Note that we don't need to integrate separately for M0 and M1. The integral calculates the appropriate expectation for individuals based on their GED status. */
			
				/*----------------------------*/
				/* Set Bounds for Integration */
				/*----------------------------*/
				
					local y1_i = y1[`i'];
					local y2_i = y2_2regime[`i'];
					local rho0 = rho0[`i'];
					local rho1 = rho1[`i'];
										
					/* Set input for determining integration bounds */
										
						local V_TruncVal = zg[`i'];       local Upper_V_TruncInd = `y1_i'; /* Note: Upper_V_TruncInd = GED corresponds to GED=1 <=> zg-v>0 <=> v<zg which is indeed an indication of upper truncation at zg. */ 
						local M_TruncVal = zf`y1_i'[`i']; local Upper_M_TruncInd = `y2_i';
					
					/* Set integration bounds for our cdf transformations -- u_v and u_m -- of their corresponding normal variables */ 
					
						if `Upper_V_TruncInd' == 1 {; local a_v = 0.0;                  local b_v = normal(`V_TruncVal'); local v_mult =  1; };
						if `Upper_V_TruncInd' == 0 {; local a_v = normal(`V_TruncVal'); local b_v = 1.0;                  local v_mult = -1; };
						if `Upper_M_TruncInd' == 1 {; local a_m = 0.0;                  local b_m = normal(`M_TruncVal'); local m_mult =  1; };
						if `Upper_M_TruncInd' == 0 {; local a_m = normal(`M_TruncVal'); local b_m = 1.0;                  local m_mult = -1; };
					
					* di "a_v = `a_v', b_v = `b_v', a_m = `a_m', b_m = `b_m'.";
					
					/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
					/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
					/*  NOTE explaining @_mult terms above: it is necessary to be creative with calling the binormal() cdf, because it only cumulates to an upper bound, and does not permit a call to calculate  
						density of the upper tail. To calculate the density above a lower bound, I make use of the symmetry of the distribution, where the area above our desired upper bound is equivalent to
						the area below the negative of that upper bound. This is equivalent to examining properties of the negative of that random variable. When using that transformation, it is necessary to
						reverse the sign of the rho, since the sign of the covariance between the two random variables is now reversed. Put more simply: 
							
								Pr(A<a, B>b, rho(A,B)) = Pr(A<a, -B<-b, rho(A,-B)) = binormal(a,-b,-rho(A,B)) 
								
									and
									
								Pr(A>a, B>b, rho(A,B)) = Pr(-A<-a, -B<-b, rho(-A,-B)) = binormal(-a,-b,rho(A,B))
								
						So, these @_mult values will be used to multiply the bound and rho in the bivariate cdf function below.   */
					/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
					/*-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/

					/* Set up Variable Transformation for Gaussian Sum. */
						/* See e.g. http://en.wikipedia.org/wiki/Gaussian_quadrature. */
				
						local t1_v = ((`b_v')-(`a_v'))/2; local t2_v = ((`b_v')+(`a_v'))/2;
						local t1_m = ((`b_m')-(`a_m'))/2; local t2_m = ((`b_m')+(`a_m'))/2;
						* di "Transformation values are `t1_v', `t2_v', `t1_m', `t2_m'.";
					

				/*---------------------*/			
				/* Run the Integration */
				/*---------------------*/

					local QSum = 0;
					local QSum_Marg = 0;
					forvalues p1 = 1(1)10 {;
					
						/* Calculate marginal expectations of each variable, i.e. E[V|V>a] and E[M|M>b] */
						local ExpectVar_i_marg = invnormal((`t1_`var'')*(`x`p1'') + (`t2_`var''));
						
						* di "Quadrature Sum Assignment: ...QSum_Marg = `QSum_Marg' + `w`p1''*( (0.5)*invnormal(`ExpectVar_i_marg') )";
						local QSum_Marg = `QSum_Marg' + `w`p1''*( (0.5)*`ExpectVar_i_marg' );
					
						forvalues p2 = 1(1)10 {;	
					
							local V_i = invnormal((`t1_v')*(`x`p1'') + (`t2_v'));
							local M_i = invnormal((`t1_m')*(`x`p2'') + (`t2_m')); /* Note that this "M" variable generically stands for either "M0" or "M1", based on what the individual's GED status is. */
							
							local binormpdf_num = exp( (-1) * (((`V_i')^2) + ((`M_i')^2) - (2*(`rho`y1_i'')*(`V_i')*(`M_i'))) / (2*(1-(`rho`y1_i'')^2)) );  /* Numerator of the bivariate normal */
							local binormpdf_den = 2*_pi*sqrt(1-(`rho`y1_i'')^2);                                                                            /* Denominator of the bivariate normal */
							local binormpdf = `binormpdf_num'/`binormpdf_den';
							
							local Summand =  (``ExpectVar'_i'*`binormpdf') / (normalden(`V_i')*normalden(`M_i'));
							local QSum = `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' );
							
							/* Misc auditing code:
								* di "V_i = `V_i', M_i = `M_i'";
								* di "binormpdf = " `binormpdf_num'/`binormpdf_den';
								* di "Summand = " (``ExpectVar'_i'*`binormpdf') / (normalden(`V_i')*normalden(`M_i'));
								* di "QSum = " `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' );
							*/
						
						};
					};
					
					local binorm_prob = binormal((`v_mult')*(`V_TruncVal'), (`m_mult')*(`M_TruncVal'), (`v_mult')*(`m_mult')*(`rho`y1_i''));
					local QSum = (`t1_v')*(`t1_m')*(`QSum')/`binorm_prob';
					qui replace E`ExpectVar'_CtrlFn = `QSum' if _n == `i';
					* di "Marginal Expectation Replacement ... qui replace E`ExpectVar'_MargExpect = `QSum_Marg' if _n == `i'"; 
					qui replace E`ExpectVar'_MargExpect = `QSum_Marg' if _n == `i';
					* pause;
				
			}; /* End of run across variables of integration */
			
		}; /* End of run across individuals */
		
		di "Integration took from `start_time' to " c(current_time);
		
		cap drop EM0_CtrlFn EM1_CtrlFn EM0_MargExpect EM1_MargExpect;
		gen EM0_CtrlFn = 0; gen EM1_CtrlFn = 0; gen EM0_MargExpect = 0; gen EM1_MargExpect = 0;
		replace EM0_CtrlFn = EM_CtrlFn if y1 == 0 & y2_2regime == 1; label var EM1_CtrlFn "Expected value of work error, conditional on GED and Work.";
		replace EM1_CtrlFn = EM_CtrlFn if y1 == 1 & y2_2regime == 1; label var EM0_CtrlFn "Expected value of work error, conditional on no GED and Work.";
		
		replace EM1_MargExpect = EM_MargExpect if y2_2regime == 1; label var EM0_MargExpect "Expected value of work error, conditional on Working.";
		replace EM0_MargExpect = EM_MargExpect if y2_2regime == 0; label var EM0_MargExpect "Expected value of work error, conditional on Not Working.";
		
		
	/*------------------------------------------------------------*/
	/* Check Quadrature Calculations with Monte Carlo Integration */
	/*------------------------------------------------------------*/

		/* Running integration for the individual in row 1,  who has GED = 0, Employment = 0, zg = -0.7001294, zf0 = -1.157879,  EV_CtrlFn = 0.4444327,  EM_CtrlFn = 0.343767. */
		/* Running integration for the individual in row 10, who has GED = 0, Employment = 1, zg = -1.31399,   zf1 = -0.02880334, EV_CtrlFn = -0.2087638, EM_CtrlFn = -0.9259795 */
		/* The variances of the random variables underlying these outcomes is assumed to be 1, and rho0 is estimated as : 0.5124192, and rho1 = 0.3725314 */
		
		/*
		#delimit;
			set more off, perm;
			cap set matsize 10000;
			
			local nDraws = 10000;
			matrix MCMat = J(`nDraws',2,.);
			
			/* Set Draw Parameters */
				/* Model is:
				u1, u2 ~ N(0,1)
				e1 = a1*u1
				e2 = a2*u1 + a3*u2
				e ~ N([0 0]', [ a1^2        .      ]
						      [ a1*a2 (a2^2 + a3^2)]) */
				local v1  = 1.0;
				local v2  = 1.0;
				local rho = 0.5124192;
				local ltrunc1 = -1.31399; /* Model is: y1 = 1 <=> zg - v > 0 <=> v < zg. For someone with y1 = 0, we want a lower bound. */
				local utrunc2 = -0.2880334;
				
				scalar a1 = `v1';
				scalar a2 = `rho'/scalar(a1);
				scalar a3 = sqrt(`v2' - scalar(a2)^2);
				
			/* Run Integration */
			
				local i = 0;
				while `i' < `nDraws' {;
					di "i = `i'";
					scalar u1 = rnormal();
					scalar u2 = rnormal();
					scalar V = scalar(a1)*scalar(u1);
					scalar M = scalar(a2)*scalar(u1) + scalar(a3)*scalar(u2);
					
					if (scalar(V) > `ltrunc1') & (scalar(M) < `utrunc2') {;
						local i = `i' + 1;
						scalar list V M;
						matrix MCMat[`i',1] = scalar(V);
						matrix MCMat[`i',2] = scalar(M);
					};
				};
				matrix MeanVec = J(1,`nDraws',1/`nDraws');
				matrix Mean = MeanVec*MCMat;
				matrix list Mean;
			
			*/
			

	/*---------------------------------------------------------------------------*/
	/* Simple Estimation Checks: Single Endogenous Variable, or Single Selection */
	/*---------------------------------------------------------------------------*/

		#delimit cr;

		/* y1 is an endogenous treatment */
		
			reg y4_endogtreat x4 y1
			cap drop y1_haz
			treatreg y4_endogtreat x4, treat(y1 = x1) hazard(y1_haz) /* First stage is right */
			reg y4_endogtreat x4 y1 y1_haz
			reg y4_endogtreat x4 y1 EV_MargExpect /* Almost identical results */
			
			/* Other procedures...
				*ivregress 2sls y4_endogtreat x4 (y1 = x1) 
				*ivregress liml y4_endogtreat x4 (y1 = x1)
				*ivregress gmm  y4_endogtreat x4 (y1 = x1)
			*/
			
			/* Recall that ...
				y1            = -1.5 + 3.0*x1 + u1 > 0
				y4_endogtreat = 0.0 + 3.0*x4 + 5.0*y1 + u4 ... where u4 is correlated with u1 (which determines y1) */
			
			
		/* y2_1regime determines observation of y4_sel_1regime */
		
			reg y4_sel_1regime x4
			cap drop y2_1regime_mills
			heckman y4_sel_1regime x1 x4, select(y2_1regime = x2) mills(y2_1regime_mills)
			reg y4_sel_1regime x4 EM_MargExpect /* Almost identical results */
		
			/* Recall that ...
				y4_sel_1regime = 0.0 + 3.0*x4 + u4 ... which is observed only if y2_1regime == 1
				y2_1regime  = ( 0.3 - 2.0*x2 + u2 > 0) ... where u2 is correlated with u4 
				... because x4 is not correlated with the error u2, this should just result in a biased estimate of the intercept term in the absence of correction.   */


	/*----------------------------------------------------------------*/		
	/* Run Outcome Estimation With Selection and Endogenous Regressor */
	/*----------------------------------------------------------------*/			

		/* Make sure that all of the control function variables only take on values for the case that they apply to. */
		#delimit;
		cap drop E*_CtrlFn_G*; 
		
		gen EV_CtrlFn_G0  = EV_CtrlFn * (y1 == 0) * (y2_2regime == 1);
		gen EV_CtrlFn_G1  = EV_CtrlFn * (y1 == 1) * (y2_2regime == 1);
		gen EM0_CtrlFn_G0 = EM_CtrlFn * (y1 == 0) * (y2_2regime == 1);
		gen EM1_CtrlFn_G1 = EM_CtrlFn * (y1 == 1) * (y2_2regime == 1);
		
		/*
			reg y3_1regime x4 y1;
			reg y3_1regime x4 y1 EV_CtrlFn;
			reg y3_1regime x4 y1           EM0_CtrlFn EM1_CtrlFn;
			reg y3_1regime x4 y1 EV_CtrlFn EM0_CtrlFn EM1_CtrlFn;
		*/
		
		reg y3_2regime x4 y1;
		reg y3_2regime x4 y1 EV_CtrlFn;
		reg y3_2regime x4 y1           EM0_CtrlFn EM1_CtrlFn;
		reg y3_2regime x4 y1 EV_CtrlFn EM0_CtrlFn EM1_CtrlFn;
		reg y3_2regime x4 y1 EV_CtrlFn_G0 EV_CtrlFn_G1  EM0_CtrlFn_G0  EM1_CtrlFn_G1;
		
		/* Recall that ... y3_1regime = 0.0 + 3.0*x4 + 5.0*y1 + u4 ... which is observed only if y2_2regime = 1 */
		
			
			
#delimit cr;
