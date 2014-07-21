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

