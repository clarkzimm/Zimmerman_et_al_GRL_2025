Code repository for 
"Slowed response of AMOC not a robust signal of collapse"
Zimmerman, C. C., Wagner, T. J. W., Maroon, E. A., & McNamara, D. E. (2025).  
Geophysical Research Letters, 51, e2024GL112415.
https://doi.org/10.1029/2024GL112415

Contact:
cczimmerman3@wisc.edu or 
till.wagner@wisc.edu

9 January 2025


All model output needed to recreate figures in the main text and S.I. provided in the folder 'outputs'. 
Scripts in folder 'model_code' used to create output files contain 'save()' functions indicating which output files were created in each script, and are additionally listed here. 
Code needed to recreate main text and S.I figures are provided in folder 'figure_code'. 
A flowchart indicating all script dependencies needed to produce each figure is provided ('Code_dependencies_flowchart.pdf').

All scripts are Matlab '.m' files and output is in Matlab '.mat' format. 

XX Dependencies are provided in the individual scripts. 

Scripts in folder 'model_code' can be used as follows:

Analytic_2box.m : This file produces the analytic equilibrium solutions of AMOC strength Q as a function of hosing H for weak or strong gyres in the 2-box model with choice of model calibration.
		Provided output: 'analytic_sol_2box_FMSB_1CO2_strong.mat'; 'analytic_sol_2box_FMSB_1CO2_weak.mat'

Analytic_3box.m : This file produces the analytic equilibrium solutions of AMOC strength Q as a function of hosing H for weak or strong gyres in the 3-box model with choice of model calibration.
		Provided output: 'analytic_sol_3box_FMSB_1CO2_strong.mat'; 'analytic_sol_3box_FMSB_1CO2_weak.mat'

Analytic_3box_rampKN.m : This file produces the analytic equilibrium solutions of AMOC strength Q as a function of North Atlantic gyre transport strength KN in the 3-box model with choice of model calibration.
		Provided output: 'analytic_sol_3box_FMSB_1CO2_rampKN.mat' 

Analytic_DeltaH.m : This file produces the analytic equilibrium hysteresis (DH) of the AMOC 3-box model in hosing (H) as a function of Atlantic gyre exchange strengths KN and KS, and boundary between monostable and bistable (KN,KS) regimes.
		Provided output: 'analytic_DH_3box_FMSB_1CO2.mat'; 'analytic_DH_3box_FMSB_2CO2.mat'; 'analytic_DH_3box_HGEM_1CO2.mat'; 'analytic_DH_3box_HGEM_2CO2.mat'; 'analytic_DH_5box_FMSB_1CO2.mat'; 'analytic_DH_3box_FMSB_1CO2_short.mat'

Det_DH_3box.m : This file produces the hysteresis width (Delta H*) in hosing (H) for deterministic runs of the 3-box model for varying gyre strengths (KN,KS).
		Provided output: 'det_DH_3box_FMSB_1CO2.mat'

IC_Kdep_3box.m : This file produces a matrix of steady state salinities (SN,ST) for initial hosing value (hmin) and varying gyre strengths (KN,KS). Does not run independently.
		Provided output: 'IC_Kdep_3box_FMSB_1xCO2_dens100.mat'

Kspace_3box_stats : This file calculates the slopes, and normalized slopes, of the autocorrelation and variance with increasing hosing of Q timeseries in the three 3-box model runs for varying gyre strengths (KN,KS).
		Provided output: 'stats_Kdep_3box_FMSB_1CO2_lag50.mat'

parameters_Xbox_YY_ZCO2.m : Files containing parameter values for the X = [3 5]-box model as calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = [1 2]xCO2. Default is 'parameters_3box_FMSB_1CO2.m'.

Qdet_Kdep_3box.m : This file computes the deterministic trajectory of AMOC strength in the 3-box model for varying gyre strengths (KN,KS). Does not run independently.
		Provided output: 'Qdet_Kdep_3box_FMSB_1xCO2_dens100_lag10.mat'; 'Qdet_Kdep_3box_FMSB_1xCO2_dens100_lag50.mat'

solve_H_3box.m : This file produces analytic equilibrium solutions of hosing (H) as a function of AMOC transport (q), and gyre transport strengths (KN,KS) for 3-box model. Does not run independently.

solve_H_5box.m : This file produces analytic equilibrium solutions of hosing (H) as a function of AMOC transport (q), and gyre transport strengths (KN,KS) for 5-box model. Does not run independently.

solve_initial_salinities.m : This file solves for the steady state salinities (ST,SN) for initial hosing value (hmin) and final hosing value (hmax) given set parameterization/gyres for 3-box model.  Does not run independently.

solve_initial_salinities_2box.m : This file solves for the steady state salinities (ST,SN) for initial hosing value (hmin) and final hosing value (hmax) given set parameterization/gyres for 2-box model. Does not run independently.

solve_initial_salinities_KN.m : This file solves for the steady state salinities (ST,SN) for initial gyre value (knmax) and final gyre value (knmin) given set parameterization/hosing for 3-box model. Does not run independently.

Stochastic_3box_DHstar.m : This file produces the hysteresis width (Delta H*) in hosing (H) for stochastic runs of the 3-box model for varying gyre strengths (KN,KS).
		Provided output: 'stoch_DH_3box_FMSB_1CO2.mat'

Stochastic_sims_2box.m : This file integrates the 2-box model for both forward and backward hosing (H), solving for AMOC strength (Q), with Guassian white noise added to the hosing, for weak or strong gyres.
		Provided output: 'stochastic_2box_FMSB_1CO2_strong.mat'; 'stochastic_2box_FMSB_1CO2_weak.mat'

Stochastic_sims_3box.m : This file integrates the 3-box model for both forward and backward hosing (H), solving for AMOC strength (Q), with Guassian white noise added to the hosing, for weak or strong gyres.
		Provided output: 'stochastic_3box_FMSB_1CO2_strong.mat'; 'stochastic_3box_FMSB_1CO2_weak.mat'

Stochastic_sims_3box_rampKN.m : This file integrates the 3-box model for both forward and backward ramping of North Atlantic gyre transport strength (KN), solving for AMOC strength (Q), with Guassian white noise added to KN
		Provided output: 'stochastic_3box_FMSB_1CO2_rampKN.mat'

An additional output: 'Qsol4_mathematica.mat' is provided in the 'outputs' folder and is called by 'Fig2_stoch_3box_stats.m', 'SI_FigS2_ratedep.m', and 'SI_FigS3_noisedep.m' to fill in the negative Q unstable analytic equilibrium solution which is not produced by MATLAB's vpasolve function. 



