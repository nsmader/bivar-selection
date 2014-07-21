/*----------------------------------------------------------------------------------*/
/* SANDBOX CODE FOR NUMERICAL INTEGRATION - BOTH MONTE CARLO AND QUADRATURE METHODS */
/*----------------------------------------------------------------------------------*/

	set matsize 10000

/*-----------------------------------------------------*/
/*-----------------------------------------------------*/
/* Monte Carlo Integration - One Dimension - Normal RV */
/*-----------------------------------------------------*/
/*-----------------------------------------------------*/

	/* Procedure:
		o Draw random variables according to the known or estimated distribution;
		o Form quantity of interest;
		o Perform many draws in order to approximate the true value of the object of interest; 
	*/
	
	/* Play-test Procedure with Known Integral */
	
		/* Numerically calculate E[e|e>1] where e~N(.25,1). */
			local mu = .25
			local v  = 1
			local s  = sqrt(`v')
			local lower = 1
			
			local z = (`lower' - `mu')/`s'

			local nDraws = 10000
			matrix MCMat = J(`nDraws',1,.)
			local p_LowerTrunc = normal(`z')
			* di `p_LowerTrunc'
			forvalues i = 1(1)`nDraws' {
				local ui = runiform()
				local ui_trunc = `p_LowerTrunc' + `ui'*(1-`p_LowerTrunc')
				local ei = invnormal(`ui_trunc')*`s' + `mu'
				matrix MCMat[`i',1] = `ei'
			}
			matrix MeanVec = J(1,`nDraws',1/`nDraws')
			matrix Mean = MeanVec*MCMat
			matrix list Mean
			
			local solution = (normalden(`z')/(1 - normal(`z')))*`s' + `mu'
			di "Analytical solution is: `solution'"
		
		
		/* Numerically calculate E[e|e<2] where e~N(1,.5). */
			local mu = 1
			local v  = .5
			local s  = sqrt(`v')
			local upper = 2
			
			local z = (`upper' - `mu')/`s'
			
			local nDraws = 10000
			matrix MCMat = J(`nDraws',1,.)
			local p_UpperTrunc = normal(`z')
			* di `p_UpperTrunc'
			forvalues i = 1(1)`nDraws' {
				local ui = runiform()
				local ui_trunc = `ui'*`p_UpperTrunc'
				local ei = invnormal(`ui_trunc')*`s' + `mu'
				matrix MCMat[`i',1] = `ei'
			}
			matrix MeanVec = J(1,`nDraws',1/`nDraws')
			matrix Mean = MeanVec*MCMat
			matrix list Mean
			
			local solution = (-normalden(`z')/normal(`z'))*`s' + `mu'
			di "Analytical solution is: `solution'"
			
		

/*-------------------------------------------------------*/
/*-------------------------------------------------------*/
/* Monte Carlo Integration - Two Dimensional - Normal RV */
/*-------------------------------------------------------*/
/*-------------------------------------------------------*/

	/* Model is:
		u1, u2 ~ N(0,1)
		e1 = a1*u1
		e2 = a2*u1 + a3*u2
		e ~ N([0 0], [ a1^2        .      ]
				     [ a1*a2 (a2^2 + a3^2)])
	*/


	/* Numerically calculate Phi12(1.0,0.5,0.2) */
		cap set matsize 10000
		
		local nDraws = 10000
		matrix MCMat = J(`nDraws',1,.)
		local v1  = 1.0
		local v2  = 1.0
		local rho = 0.2
		local trunc1 = 1.0
		local trunc2 = 0.5
		
		local p = binormal(`trunc1', `trunc2', `rho')
		di `p'
		
		local a1 = `v1'
		local a2 = `rho'/`a1'
		local a3 = sqrt(`v2' - `a2'^2)
		
		forvalues i = 1(1)`nDraws' {
			local u1 = rnormal()
			local u2 = rnormal()
			local e1 = `a1'*`u1'
			local e2 = `a2'*`u1' + `a3'*`u2' /* NSM: Note... This code is still an unfinished application of the GHK sampler. I am reading notes from Kenneth Train's book on this. */
			
			local I = (`e1' < `trunc1') & (`e2' < `trunc2')
			matrix MCMat[`i',1] = `I'
		}
		matrix MeanVec = J(1,`nDraws',1/`nDraws')
		matrix Mean = MeanVec*MCMat
		matrix list Mean
		
		local solution = binormal(1.0, 0.5, 0.2)
		di "Analytical solution is: `solution'"
	
	
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/* Run Gaussian Quadrature - One Dimension - Integral of Polynomial */
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/


	/* Quadrature weights are available from: http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php */

	#delimit;
	local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
	local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 

	local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
	local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;
			
	/* Calculate the integral for f(x) = 3.0*x^2 - 2.0*x + 1.0 from -1 to 3. The solution equals [x^3 - x^2 + x + k] evaluated at 3 and -1, yielding [27 - 9 + 3 + k] - [(-1) - (+1) + (-1) + k] = 24. 
		The calculation for this should be exact. (Some other calculations: IntFn(3) = 21, IntFn(2) = 6, IntFn(1) = 1, IntFn(0) = 0, IntFn(-1) = -1.)
		
		To perform this integral, we use xt which has bounds -1 to 1, by plugging values of x into the equation: xt = t1*x + t2 
	*/  
		
		local a = -1; local b = 3; 
		local t1 = ((`b')-(`a'))/2;  local t2 = ((`b')+(`a'))/2;
		di "a = `a', b = `b', t1 = `t1', t2 = `t2'";
			
		local QSum = 0;
		forvalues p = 1(1)10 {;
		
			local xt   = `t1'*`x`p'' + `t2';
			
			local QSum = `QSum' + (`w`p'')*((3.0*(`xt')^2 - 2.0*(`xt') + 1.0)*(`t1'));
			
		};
		di `QSum';

di normalden(2,1,.5)
di normalden((2-1)/.5)/.5


/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/* Run Gaussian Quadrature - One Dimension - Integral of Normal RV */
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

	#delimit;
	local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
	local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 

	local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
	local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;
	
	/* Calculate the expected value E[V|V>lower] where V~N(mu,v) */
	 
		/* For some notes explaining the implementation of the quadrature procedure, see the "Notes on Gaussian Quadrature.lyx" file in my "writeup" folder. */  
		
		local mu = 2;
		local v  = 3;
		local s  = sqrt(`v');
		local lower = 1.35;
		
		local z = (`lower' - `mu')/`s'; /* This is the lower bound of the standard normal transformation of our random variable */
		di `z';
		
		local a = normal(`z'); local b = 1; 
		local t1 = ((`b')-(`a'))/2;
		local t2 = ((`b')+(`a'))/2;
		di "b = `b', a = `a', t1 = `t1', t2 = `t2'";
			
		local QSum = 0;
		forvalues p = 1(1)10 {;
		
			local xt   = `t1'*`x`p'' + `t2';
			local QSum = `QSum' + `w`p''*( (0.5)*invnormal(`xt') );
			
		};
		local QSum = `s'*`QSum' + `mu';
		di `QSum';
		
		local solution = (normalden(`z')/(1 - normal(`z')) )*`s' + `mu';
		di "Analytical solution is: `solution'";


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/* Run Gaussian Quadrature - Two Dimensions - Integral of Polynomial */
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

	#delimit;
	local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
	local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 

	local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
	local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;
	
	/* Calculate the double integral for f(x,y) = 4.0*xy for given bounds: */
	
		/* Test 1: Integral of (x: 2-3, y: 1-3) is... Integral(2.0*x^2*y*dy...evaluated at 2,3) = Integral(10.0*y*dy) = 5.0*y^2...evaluated at 1,3 = 5.0*(9) - 5.0(1) = 40.0 */
		
		/* Test 2: Integral of (x: 0-3, y: 0-1) is... Integral(2.0*x^2*y*dy...evaluated at 0,3) = Integral(18.0*y*dy) = 9.0*y^2...evaluated at 0,1 = 9.0*(1) - 9.0(0) = 9.0 */
			
		local ax = 0; local bx = 3;
		local ay = 0; local by = 1;
		
		local t1x = ((`bx')-(`ax'))/2; local t2x = ((`bx')+(`ax'))/2;
		local t1y = ((`by')-(`ay'))/2; local t2y = ((`by')+(`ay'))/2;
		di "Transformation values are `t1x', `t2x', `t1y', `t2y'.";
			
		local QSum = 0;
		forvalues p1 = 1(1)10 {;
			forvalues p2 = 1(1)10 {;	
		
				local xt1 = (`t1x')*(`x`p1'') + (`t2x');
				local xt2 = (`t1y')*(`x`p2'') + (`t2y');
				
				* di "QSum cmd is ... local QSum = `QSum' + (`w`p1'')*(`w`p2'')*( 4.0*`xt1'*`xt2' )"; 
				local QSum = `QSum' + (`w`p1'')*(`w`p2'')*( 4.0*(`xt1')*(`xt2')*(`t1x')*(`t1y'));
			
			};
		};
		local QSum = `QSum';
		di `QSum';



/*-----------------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------*/
/* Run Gaussian Quadrature - Two Dimensions - Expected Value of Normal RV, Conditioning on BVN Density */
/*-----------------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------------*/


	/*--------------------------------------------------------------------*/
	/* Calculate the expected value E[V|V>a, M>b] where V,M ~ N(0,[1 . ]
															      [r 1 ])
	/*--------------------------------------------------------------------*/
		/* See the notes on the integrand which in the file "../writeups/Notes on Gaussian Quadrature.lyx"  */
	
	
		#delimit;
		local x1 -0.97390653; local x2 -0.86506337; local x3 -0.67940957; local x4 -0.43339539; local x5 -0.14887434; 
		local x10 0.97390653; local x9  0.86506337; local x8  0.67940957; local x7  0.43339539; local x6  0.14887434; 
	
		local w1  0.06667134; local w2  0.14945135; local w3  0.21908636; local w4  0.26926672; local w5  0.29552422;
		local w10 0.06667134; local w9  0.14945135; local w8  0.21908636; local w7  0.26926672; local w6  0.29552422;
					
		local V_TruncVal = 0.2; local Lower_V_Trunc = 0; /* This is chosen to correspond to the outcome of interest, e.g. Lower_V_Trunc = (GED==1), since Pr(GED=1|z) = Pr(zg+v>0) = Pr(v>-zg) is a lower truncation. */ 
		local M_TruncVal = 0.5; local Lower_M_Trunc = 1;
		local rho = 0.3;
		
		if `Lower_V_Trunc' == 1 {; local a_v = normal(`V_TruncVal'); local b_v = 1.0;                  local biv_mult_v = -1; };
		if `Lower_V_Trunc' == 0 {; local a_v = 0.0;                  local b_v = normal(`V_TruncVal'); local biv_mult_v =  1; };
		if `Lower_M_Trunc' == 1 {; local a_m = normal(`M_TruncVal'); local b_m = 1.0;                  local biv_mult_m = -1; };
		if `Lower_M_Trunc' == 0 {; local a_m = 0.0;                  local b_m = normal(`M_TruncVal'); local biv_mult_m =  1; };
		
		di "a_v = `a_v', b_v = `b_v', a_m = `a_m', b_m = `b_m'.";
		
		/* NOTE explaining biv_mult_@ terms above: it is necessary to be creative with calling the binormal() cdf, because it only cumulates to an upper bound, and does not permit a call to calculate  
			density of the upper tail. To calculate the density above a lower bound, I make use of the symmetry of the distribution, where the area above our desired upper bound is equivalent to
			the area below the negative of that upper bound. This is equivalent to examining properties of the negative of that random variable. When using that transformation, it is necessary to
			reverse the sign of the rho, since the sign of the covariance between the two random variables is now reversed. Put more simply: 
				
					Pr(A<a, B>b, rho(A,B)) = Pr(A<a, -B<-b, rho(A,-B)) = binormal(a,-b,-rho(A,B)) 
					
						and
						
					Pr(A>a, B>b, rho(A,B)) = Pr(-A<-a, -B<-b, rho(-A,-B)) = binormal(-a,-b,rho(A,B))
					
			So, these biv_mult_@ values will be used to multiply the bound and rho in the bivariate cdf function below.	
		*/  
		
		local t1_v = ((`b_v')-(`a_v'))/2; local t2_v = ((`b_v')+(`a_v'))/2;
		local t1_m = ((`b_m')-(`a_m'))/2; local t2_m = ((`b_m')+(`a_m'))/2;
		di "Transformation values are `t1_v', `t2_v', `t1_m', `t2_m'.";
			
		local QSum = 0;
		forvalues p1 = 1(1)10 {;
			forvalues p2 = 1(1)10 {;
		
				local V = invnormal((`t1_v')*(`x`p1'') + (`t2_v')); 
				local M = invnormal((`t1_m')*(`x`p2'') + (`t2_m'));
				
				local binormpdf_num = exp( (-1) * (((`V')^2) + ((`M')^2) - (2*(`rho')*(`V')*(`M'))) / (2*(1-(`rho')^2)) ); /* Numerator of the bivariate normal expression */
				local binormpdf_den = 2*_pi*sqrt(1-(`rho')^2) ;                                                            /* Denominator of the bivariate normal expression */
				local binormpdf = `binormpdf_num'/`binormpdf_den';				
				
				local Summand =  ((`V')*(`binormpdf')) / (normalden(`V')*normalden(`M'));   /* The denominator here comes from the change of variables from integration of V and M to percentiles */
				
				local QSum = `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' );
				
				/* Misc auditing code:
					* di "V = `V', M = `M'.";
					* di "binormpdf_num calc is: exp( (-1) * (((`V')^2) + ((`M')^2) - (2*(`rho')*(`V')*(`M'))) / (2*(1-(`rho')^2)) )";
					* di "binormpdf_num = `binormpdf_num', binormpdf_den = `binormpdf_den', binormpdf = `binormpdf'.";
					* di "SUM TERM ... local Summand =  ((`V')*(`binormpdf')) / (normalden(`V')*normalden(`M'))";
					* di "QSum = `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' ) = " `QSum' + (`w`p1'')*(`w`p2'')*( `Summand' );
					* di "Summand =  ((`V')*`binormpdf') / (normalden(`V')*normalden(`M')) = `Summand'";
					* di "Summand = `Summand', QSum = `QSum'.";
					* di "QSum =  `QSum'";
				*/
			
			};
			
		};
		
		di "binormal(`biv_mult_v'*`V_TruncVal', `biv_mult_m'*`M_TruncVal', (`biv_mult_v')*(`biv_mult_m')*`rho') = " 
			binormal(`biv_mult_v'*`V_TruncVal', `biv_mult_m'*`M_TruncVal', (`biv_mult_v')*(`biv_mult_m')*`rho');
			
		local QSum = (`t1_v')*(`t1_m')*`QSum'/binormal((`biv_mult_v')*(`V_TruncVal'), (`biv_mult_m')*(`M_TruncVal'), (`biv_mult_v')*(`biv_mult_m')*`rho');
		di `QSum';
		
		
		/*---------------------------------------------------------*/
		/* Check the above solution with a MC Integration Exercise */
		/*---------------------------------------------------------*/
		
			/* Numerically calculate E[V|V>1.0,M>0.5] where V,M~MVN(0.0,v1=1.0,v2=1.0,rho=0.5) */
			#delimit;
			set more off, perm;
			cap set matsize 10000;
			
			local nDraws = 10000;
			matrix MCMat = J(`nDraws',1,.);
			
			/* Set Draw Parameters */
				/* Model is:
				u1, u2 ~ N(0,1)
				e1 = a1*u1
				e2 = a2*u1 + a3*u2
				e ~ N([0 0]', [ a1^2        .      ]
						      [ a1*a2 (a2^2 + a3^2)]) */
				local v1  = 1.0;
				local v2  = 1.0;
				local rho = 0.3;
				local utrunc1 = 0.2;
				local ltrunc2 = 0.5;
				
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
					
					if (scalar(V) < `utrunc1') & (scalar(M) > `ltrunc2') {;
						local i = `i' + 1;
						scalar list V;
						matrix MCMat[`i',1] = scalar(V);
					};
				};
				matrix MeanVec = J(1,`nDraws',1/`nDraws');
				matrix Mean = MeanVec*MCMat;
				matrix list Mean;
															       
