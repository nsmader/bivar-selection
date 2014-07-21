/*------------------------------------*/
/* SANDBOX CODE FOR RUNNING STATA MLE */
/*------------------------------------*/

/* TO DO  
	o Try estimation using d1 method (which uses a specified gradient in its search function)
	o Determine a way to report transformed rho estimates at the end of the procedure
	o (For later) write a joint estimation of the probability of GED certification and work conditional on GED status 
*/


/*-----------------------*/
/* Generate Toy Data Set */
/*-----------------------*/

	clear all
	set mem 1g
	set more off

	/*--------------------*/
	/* Generate Variables */
	/*--------------------*/
	
		set obs 10000
		set seed 1092011
	
		gen x1 = runiform()
		gen x2 = runiform()
		gen x3 = runiform()
		gen x4 = runiform()
	
	/*-----------------*/
	/* Generate Errors */
	/*-----------------*/
	
		/* Can find notes for this at http://www.stata.com/support/faqs/stat/mvnorm.htmlhttp://www.stata.com/support/faqs/stat/mvnorm.html */
	
		gen c1 = invnorm(runiform())
		gen c2 = invnorm(runiform())
		gen c3 = invnorm(runiform())
		gen c4 = invnorm(runiform())
		
		matrix Sig_u = (1.0, 0.5, 0.0, 0.9 \ ///
		                0.5, 1.0, 0.0, 0.3 \ ///
		                0.0, 0.0, 1.0, 0.3 \ ///		                
		                0.9, 0.3, 0.3, 1.0 )
		             
		matrix Chol_Sig_u = cholesky(Sig_u)
		
		local rows = rowsof(Chol_Sig_u)
		
		foreach r of numlist 1(1)`rows' {
			matrix CSu`r' = Chol_Sig_u[`r',1...]
			matrix score u`r' = CSu`r' /* This works because the default matrix names are c# and that's what we've used as our variable names for orthogonal errors */
		} 
 
 	/*-------------------*/
	/* Generate Outcomes */
	/*-------------------*/

		gen     y1  = (-1.5 + 3.0*x1 + u1 > 0)             /* This is whether an individual gets a GED */
		gen     y2  = ( 0.3 - 2.0*x2 + u2 > 0) if y1 == 0  /* This is whether an individual works or not, conditional on not having a GED */
		replace y2  = ( 0.4 - 1.0*x3 + u3 > 0) if y1 == 1  /* This is whether an individual works or not, conditional on     having a GED */
		
		gen     y2b = ( 0.3 - 2.0*x2 + u2 > 0)  /* This a non-regime-specific outcome for use in testing bivariate probit code */ 
		
		gen     y3 = .
		replace y3 = 0.0 + 1.0*x1 + 3.0*x4 + 5.0*y1 + u4 if y2 == 1
		

/*---------------------------*/		 
/*---------------------------*/
/* #1. RUN PROBIT ESTIMATION */
/*---------------------------*/
/*---------------------------*/

	/*----------------*/
	/* Define Program */
	/*----------------*/

		program drop _all
		
		program MyProbit
			args lf zg
			qui replace `lf' = ln(normal(`zg'))*$ML_y1 + ln(normal(-`zg'))*(1-$ML_y1)
		end

	/*-----------------*/
 	/* Call Estimation */
 	/*-----------------*/
 	
 		ml model lf MyProbit (y1 = x1)
 		ml maximize
 		
 		probit y1 x1
 		
 		/* Recall: y1  = (-1.5 + 3.0*x1 + u1 > 0)  */
 		
 		/* News: both get the same result. Good stuff. */
 		/* Note: all objects are accessible through e() objects */
 	

/*-------------------------------------*/ 	
/*-------------------------------------*/
/* #2. RUN BIVARIATE PROBIT ESTIMATION */
/*-------------------------------------*/
/*-------------------------------------*/

	/*----------------*/
	/* Define Program */
	/*----------------*/

		cap program drop MyBiProb
		program MyBiProb
		
			args lf zg1 zg2 logitrho  /* Note: logitrho equals logit((rho+1)/2) so that we represent the parameter on the real line */  
			
			tempvar rho Phi1 Phi2 Phi12 p11 p10 p01 p00
			
			qui gen double `rho' = invlogit(`logitrho')*2 - 1
			
			qui gen double `Phi1'  = normal(`zg1')                 /* = p11 + p10 */
			qui gen double `Phi2'  = normal(`zg2')                 /* = p11 + p01 */
			qui gen double `Phi12' = binormal(`zg1', `zg2', `rho') /* = p11 */
			
			qui gen double `p11' = `Phi12'
			qui gen double `p10' = `Phi1' - `Phi12'
			qui gen double `p01' = `Phi2' - `Phi12'
			qui gen double `p00' = 1 - `p11' - `p10' - `p01'
			
			qui replace `lf' = (  $ML_y1) * (  $ML_y2) * ln(`p11') ///
			                 + (1-$ML_y1) * (  $ML_y2) * ln(`p01') ///
			                 + (  $ML_y1) * (1-$ML_y2) * ln(`p10') ///
			                 + (1-$ML_y1) * (1-$ML_y2) * ln(`p00')
			
		end

	/*-----------------*/
	/* Call Estimation */
	/*-----------------*/
	
		/*----------------------*/
		/* Initialize estimates */
		/*----------------------*/
		
			ml model lf MyProbit (y1 = x1)
			ml maximize
			matrix b1 = e(b)
			
			ml model lf MyProbit (y2b = x2)
			ml maximize
			matrix b2 = e(b)
			
			corr y1 y2b
			local logitrho = logit((r(rho)+1)/2)
			
			matrix b0 = b1, b2, scalar(`logitrho')
		
		
		/*------------------------*/
		/* Run and Compare Models */
		/*------------------------*/
		
			ml model lf MyBiProb (y1 = x1) (y2b = x2) (), vce(robust)
			ml init b0, copy
			ml maximize
			
			local rho = invlogit(_b[eq3:_cons])*2-1
			di "rho = `rho'"			
			
			biprobit (y1 = x1) (y2b = x2)
			
			/* Recall: 
				y1  = (-1.5 + 3.0*x1 + u1 > 0)  
				y2b = ( 0.3 - 2.0*x2 + u2 > 0)
				rho(y1,y2b) = 0.5
			*/
		

		
/*---------------------------------------------*/
/* Run Conditional Probit Estimation - 2 stage */
/*---------------------------------------------*/

	/*----------------*/
	/* Define Program */
	/*----------------*/
	
		cap program drop CondProb_2stage
		program CondProb_2stage
		
			args lf zf0 zf1 logitrho0 logitrho1  /* Note: the logitrho terms equal logit((rho+1)/2) so that we represent the parameter on the real line */
				
			/* Note: In this function, choice 1 represents the decision to GED certify. This is not modeled here, but rather is controlled by a previous stage's estimation (thus the reason why this
				is a "2 stage" procedure and the predicted index value "zg" is stored in a data variable). Choice 2 is the choice to work or not, conditional on choice 1. The parameters governing 
				choice 2 conditional on GED=0 and GED=1 are predicted here, and are represented in the indices "zf0" and "zf1" (short for "z phi sub 0" and "z phi sub 1"). */  
			
			tempvar rho0 rho1 p00 p01 p10 p11
			
			qui gen double `rho0' = invlogit(`logitrho0')*2 - 1
			qui gen double `rho1' = invlogit(`logitrho1')*2 - 1
			
			/* Note: it is necessary to be creative with calling the binormal() cdf, because it only cumulates density below an upper bound (and does not permit a call to calculate density of the upper tail).  
				I am currently making use of the symmetry of the normal distribution, where the area above our desired upper bound is equivalent to the area below the negative of that upper bound.
				This is effectively equivalent to examining properties of the negative of that random variable. When using that transformation, it is necessary to reverse the sign of the rho,
				since the sign of the covariance between the two random variables is now reversed. Put more simply: 
				
					Pr(A<a, B>b, rho(A,B)) = Pr(A<a, -B<-b, rho(A,-B)) = binormal(a,-b,-rho(A,B)) 
					
						and
						
					Pr(A>a, B>b, rho(A,B)) = Pr(-A<-a, -B<-b, rho(-A,-B)) = binormal(-a,-b,rho(A,B))	
				
			*/
			
			/* In our model (following the terms and structure of Professor Heckman's notes on the dual selection procedure):
					
					p11 = Pr(GED=1, Empl=1|Z) = Pr(zg - V > 0, zf - M > 0) = Pr(V < zg, M < zf) = Phi_2(zg, zf(GED), rho(V,M)).
					
				The other probabilities are suitably modified according to the logic in the above comment. */
			
			qui gen double `p00'   = binormal(-$PredZG, -`zf0',  `rho0')
			qui gen double `p01'   = binormal(-$PredZG,  `zf0', -`rho0')
			qui gen double `p10'   = binormal( $PredZG, -`zf1', -`rho1')
			qui gen double `p11'   = binormal( $PredZG,  `zf1',  `rho1')
			
			/* Note that "$Ind_y1" is an indicator value of whether someone has y1=1. This is set by the declaration of a global variable to reference a variable in the data. */
			/* The "$ML_y1" and "$ML_y2" terms are automatically generated by the MLE function to be the variables declared as the dependent variables in the 1st and 2nd equations fed to the MLE. 
				In this case, these are the values of y2 (under regime GED = 0) and y2 (under regime GED = 1). The fact that they are labeled as "...y1" and "...y2" by default is unfortunate coincidence. */  
			
			qui replace `lf' = (1-$Ind_y1) * (1-$ML_y1) * ln(`p00') ///
			                 + (1-$Ind_y1) * (  $ML_y1) * ln(`p01') ///
			                 + (  $Ind_y1) * (1-$ML_y2) * ln(`p10') ///
			                 + (  $Ind_y1) * (  $ML_y2) * ln(`p11')
		end


	/*----------------------*/
	/* Initialize estimates */
	/*----------------------*/
		
		ml model lf MyProbit (y1 = x1), robust
		qui ml maximize
		cap drop zg
		predict zg, xb
		global PredZG zg
		global Ind_y1 y1
		
		ml model lf MyProbit (y2 = x2) if y1 == 0, robust
		qui ml maximize
		matrix b2_0 = e(b)
		
		ml model lf MyProbit (y2 = x3) if y1 == 1, robust
		qui ml maximize
		matrix b2_1 = e(b)
		
		corr y1 y2
		local logitrho = logit((r(rho)+1)/2)
		
		matrix b2_init = b2_0, b2_1, scalar(`logitrho'), scalar(`logitrho') 
		
	/*----------------------------------------------------*/
	/* Run and Compare Models for Dual Hiring Regime Case */
	/*----------------------------------------------------*/
	
		ml model lf CondProb_2stage (h0: y2 = x2) (h1:y2 = x3) /lrho0 /lrho1, robust  /* Call the estimation procedure */
		ml init b2_init, copy                                                         /* Specify the values for initialization */
		ml maximize                                                                   /* Set estimation to run */
		local rho0 = invlogit(_b[lrho0:_cons])*2-1                                    /* Transform and report rho estimates */
		local rho1 = invlogit(_b[lrho1:_cons])*2-1
		di "rho0 = `rho0', rho1 = `rho1'"

		/* Recall:
			y2  = ( 0.3 - 2.0*x2 + u2 > 0) if y1 == 0	
			y2  = ( 0.4 - 1.0*x3 + u3 > 0) if y1 == 1
			rho(y1,y2(y1=0)) = 0.5, rho(y1,y2(y1=1)) = 0.0 
		*/ 

 	
 	
/*------------------------------------------------------*/
/* Run Conditional Probit Estimation - Joint Likelihood */
/*------------------------------------------------------*/ 	
 	

/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/* Final stage to check the control function method to get consistent estimates in prediction of y3 are performed in the "Sandbox Code for Full Control Function Procedure.do" file */
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/


