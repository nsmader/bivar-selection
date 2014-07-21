*change directory here
cd "C:\First Year 2010-11\Summer\Heckman"

capture log close
*log using output,text replace
*Created by: Chanont (Big) Banternghansa
*Date Created: 7/27/11
*Late Modified: 8/03/11
*Description: Create a basic regime switching model of the GED


clear all
program drop _all
program switch_lf

  args lnf mu1 mu2 mu3 sig1 sig2	

  tempvar u1 u2 eta l_j
  qui {
    gen double `u1'=$ML_y1-`mu1'
    gen double `u2'=$ML_y2-`mu2'
    gen double `eta'=(`mu3'/1)

    scalar const1=0.5*ln(2*_pi*`sig1')
    scalar const2=0.5*ln(2*_pi*`sig2')

    gen double `l_j' = -const2-0.5*(`u2'^2)/`sig2'                    if $panel==1 			        /*(wage)log contribution for ppl who switched*/    
    replace  `l_j' = ln(normal( `eta'))-const2-0.5*(`u2'^2)/`sig2'    if $ML_y3 == 1 			    /*(wage and dec to switch) log contribution for ppl who JUST switched*/
    replace `l_j' = ln(normal(-`eta'))-const1-0.5*(`u1'^2)/`sig1'     if $ML_y3 == 0 & $panel==0 	/*(wage) log contribution for ppl who havent switched*/
    
    quietly replace `lnf'=`l_j'
  }
   
end


clear
eststo clear
set obs 1000
set seed 5555

*simulate model
gen x1 = invnormal(uniform())
gen x2 = invnormal(uniform())
generate z= invnorm(uniform()) 				      		/*gen z    */
generate v= invnorm(uniform()) 				      		/*gen v    */
generate d=(3*z>v)										/*gen dummy that holds GED at the time */
gen y=1+3*x1+2.5*x2+2*invnorm(uniform())				/*gen pre regime wage*/
replace y=2.5+5*x1+7*x2+invnorm(uniform()) if d==1		/*gen post regime wage*/

*generate fake people
gen grp=.
forvalues i=1/100 {
qui replace grp=`i' if _n<=(`i'*10) & _n>(`i'-1)*10
}
sort grp d
by grp: gen year=1990+_n
by grp: gen thres=(d[_n]>d[_n-1])					/*indicator for ppl who just obtained GED*/
global panel = "d"
order year grp d thres

*Perform estimation
ml model lf switch_lf ///
 (Pre_Regime:y=x1 x2) ///				
 (Post_Regime:y=x1 x2) ///
 (Switcher:thres=z,nocons) ///
 (Pre_Std:) ///
 (Post_Std:)
ml maximize
