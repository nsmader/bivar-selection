/*----------------*/
/* PANEL ANALYSIS */
/*----------------*/

#delimit;
clear all;
cap set mem 5g;
cap set maxvar 30000;
cap set matsize 100000;
cap set more off, perm;

local InputDataPath = "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Intermediate Data";
local SavedDataPath = "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Experience and Panel Estimation/data";
local MatricesPath  = "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Experience and Panel Estimation/output/matrices";
local TablesPath    = "/mnt/ide0/share/klmshare/GED/GED Book/Data Sets/NLSY79/NLSY79 Replication and Panel Estimates/Experience and Panel Estimation/output/tables";

use "`InputDataPath'/NLSY79 - Added Var Construction - Has Factors.dta", clear;

#delimit;
gsort id - year;
	foreach v of varlist unemp1-unemp4 lowage1-lowage5 {;
		replace `v' = `v'[_n-1] if `v'[_n] == . & id[_n] == id[_n-1];
	};


#delimit;
local ExpType Pot;

local start = c(current_time);
local Restr female == 1 & CrossSectSample == 1 & (dropout_now | ged_now);
cap drop idLastObs;
bys id: gen idLastObs = (_n == _N);
qui count if idLastObs == 1;
local nClusters = r(N);

matrix ProbMat = J(100, , .);
matrix OutMat  = J(100, , .);
forvalues r = 1(1)100 {;
	di `r'; di c(current_time);
	preserve;
		bsample `nClusters', cluster(id);
		qui probit ged_now afqt_pre_std_inA* famstatus_broken_inA* dregion*_inA* black_inA* hisp_inA* yearind_D* A2024 A2529 A3034 A3539 married_now_inA* SpouseWage_inA* num_children_in_hh_inA* baby_in_hh_inA* toddler_in_hh_inA* if `Restr', robust;
		cap drop zg zg_p zg_c invmills; predict zg, xb; gen zg_p = normalden(zg); gen zg_c = normal(zg); gen invmills = zg_p/zg_c;
		qui reg LvlW afqt_pre_std_inA* res_city_inA* mhgc_inA* urban_age14_inA* south_age14_inA* famstatus_broken_inA* dregion*_inA* black_inA* hisp_inA* yearind_D* A2024 A2529 A3034 A3539 invmills
			`ExpType'ExpSpl* `ExpType'ExpPostGEDSpl* if `Restr', vce(cluster id);
	restore;	
};
local end = c(current_time);
di "Start time was `start', End time was `end'."



/* Prototype code for calculating mean and var-covar matrix from the bootstrapped estimation */

		/* Attempted code post-estimation */
		matrix ones_vec = J(1,`nBoot',1/`nBoot');
		matrix b_mean = ones_vec*bBootMat; /* Contains a row vector of averaged estimates */
		matrix b_meanmat = b_mean;
		forvalues r = 1(1)`nBoot' {; matrix b_meanmat = b_meanmat \ b_meanmat; };
		matrix bBootMat_centered = bBootMat - `b_meanmat';
		matrix bVarCovar = bBootMat_centered'*bBootMat_centered;


		/* Test code for working off of data set */
		#delimit;
		keep if male == 1 & age == 30;
		local SampleN = _N; di `SampleN';
		matrix bBootMat = J(`SampleN',2,.);
		forvalues i = 1(1) `SampleN' {;
			matrix bBootMat[`i',1] = LvlW[`i'];
			matrix bBootMat[`i',2] = afqt_pre_std[`i'];
		};
		
		matrix ones_vec = J(1,`SampleN',1/`SampleN');
		matrix b_mean = ones_vec*bBootMat; /* Contains a row vector of averaged estimates */
		matrix b_meanmat = b_mean;
		forvalues r = 1(1)`SampleN' {; matrix b_meanmat = b_meanmat \ b_mean; };
		matrix bBootMat_centered = bBootMat - b_meanmat;
		matrix bVarCovar = bBootMat_centered'*bBootMat_centered;
