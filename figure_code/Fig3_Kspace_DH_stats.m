% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Creates Figure 3 from Zimmerman et al (2024). Plots the width of the
% bistable stolution (DeltaH), width of stochastic hysteresis (DeltaH*), 
% and the slope of variance and aurocorrelation as functions of gyre 
% transport strengths (KN,KS). Creates Figure S7 from the S.I. plotting the 
% normalized slopes of variance and ac. In manuscript, all for 3-box model
% calibrated to FAMOUS_B 1xCO2. Here, can adjust to other model
% formulations and parameterizations.
%
% Inputs:
%   - analytic_DH_Xbox_YY_ZCO2.mat (analytic equilibrium hysteresis (DH) in 
%       hosing (H) for varying gyre strengths (KN,KS), and boundary between 
%       monostable and bistable (KN,KS) regimes for X box model, with YY 
%       ZxCO2 parameterization)
%   - stoch_DH_3box_YY_ZCO2.mat (hysteresis width (Delta H*) in hosing (H) 
%       for stochastic runs of the 3-box model for varying gyre strengths 
%       (KN,KS), for YY ZxCO2 parameterization)
%   - stats_Kdep_3box_YY_ZCO2_lag#.mat (rise in statistical indicators
%       field as functions of gyre strengths (KN,KS) for YY ZxCO2 
%       parameterization, calculated at lag # (defined as 'int' below))
%
% Dependancies:
%   - Analytic_DeltaH.m (solves for analytic equilibrium hysteresis (DH) in 
%       hosing (H) for varying gyre strengths (KN,KS), and boundary between 
%       monostable and bistable (KN,KS) regimes if output file for
%       model/param choice does not exist yet)
%   - Stochastic_3box_DHstar.m (solves for hysteresis width (DH*) in hosing
%       (H) for stochastic simulations of the 3-box model if output file
%       for specified parameterization does not exist yet)
%   - Kspace_3box_stats.m (solves for the slope of variance and
%       autocorrelation with increasing hosing for varying gyre strengths
%       (KN,KS) if output file for specified parameterization does not 
%       exist yet)
%          
% Output:
%   - Figure 3a,b,c and d from Zimmerman et al (2024)
%   - Figure S7 from Zimmerman at al (2024)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
%close all; clear all;
%% Choose Model formulation and calibration 
% for Figure 3 in Zimmerman et al (2024) set box = 3, model = 'FMSB', CO2=1
% set model formulation (3 -> 3-box; 5 -> 5-box)
box = 3;

%choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
model = 'FMSB';
%choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
CO2 = 1;
%% Delta H

% first load in or solve for the width of the unstable analytic solution 
% associated with the model calibration of choice 
% 'analytic_DH_3box_FMSB_1CO2.mat' available in Zenodo in 'Output' folder
analytic_name = sprintf('analytic_DH_%dbox_%s_%dCO2.mat',box,model,CO2);
if exist(analytic_name,'file')==0
    Analytic_DeltaH %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_name)

%plot
figure(36); clf
levs = 0:0.025:.4; %define contour levels

cmap = colormap((brewermap(length(levs)-1,'YlGnBu'))); %define colormap, NOTE: need brewermap package to get this cmap
contourf(KNvec,KSvec,DH',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
hold on

cc = colorbar;
plot(KNvec,DHbound,'k--','linewidth',2,'HandleVisibility','off') %add boundary between monostable and bistable regimes
clim([0 0.4])

xlim([2 50])
ylim([2 50])

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH (Sv)')
%%
%add markers for scenario 1 & 2 gyre values
figure(36)
hold on
scatter(5.456,5.447,100,'m','hexagram','filled') %weak gyre
scatter(27,27,100,'c','hexagram','filled') %strong gyre
%% Delta H*
% next, load in or solve for the stochastic hysteresis associated with the 
% model calibration of choice 'stoch_DH_3box_FMSB_1CO2.mat' available in 
% Zenodo in 'Output' folder
stoch_name = sprintf('stoch_DH_3box_%s_%dCO2.mat',model,CO2); %Note: only 3-box formulation available here
if exist(stoch_name,'file')==0
    Stochastic_3box_DHstar %If output does not exist, takes ~15 hours to run on a standard laptop Nov'24
end
load(stoch_name)

%plot
figure(13);clf
%define colormap levels and colors
levs = 0:0.025:.4;
cmap = colormap((brewermap(length(levs),'YlGnBu')));

contourf(KNramp,KSramp,DH',levs,'EdgeAlpha',.25); colorbar;
cc = colorbar; 

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH* (Sv)')
clim([0 0.4])
set(cc,'ticklabels',{levs(1:2:end),0.4},'ticks',levs(1:2:end))

%Note: need to load analytic file first to plot DHbound
hold on
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4]) %add boundary between monostable and bistable regimes
%% Variance and autocorrelation
%choose your lag for the autocorrelation
% for Figure 3c,d in Zimmerman et al (2024) set int = 50
int = 10;
% load in or solve for the slopes of statistical indicators associated with 
% the model calibration and lag of choice 
% 'stats_Kdep_3box_FMSB_1CO2_lag50.mat' available in Zenodo in 'Output' 
% folder
stat_name = sprintf('stats_Kdep_3box_%s_%dCO2_lag%d.mat',model,CO2,int); %Note: only 3-box formulation available here
if exist(stat_name,'file')==0
    Kspace_3box_stats %If output does not exist, takes ~14 hours to run on a standard laptop Nov'24 for lag 50
end
load(stat_name)
%% true slopes
% for Figure 3 in Zimmerman et al (2024)

%plot
figure(9);clf
%define colormap levels and colors
levs = 0:0.01:1;
cmap = colormap(flipud((brewermap(length(levs)-1,'PuBuGn')))); %requires package 'brewermap'
%plot
contourf(KNramp,KSramp,var_slope',levs,'EdgeAlpha',0);
clim([0.15 .6])
hold on
cc = colorbar;
hold on
%add DH=0 boundary
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4]) 

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
xlim([2 50])
ylim([2 50])
ylabel(cc,'$\frac{\delta\sigma^\mathsf{2}}{\delta\overline{\mathsf{H}}} (\mathsf{Sv}^{\mathsf{-1}})$','Interpreter','latex','FontSize',12)

figure(18);clf
%set colormap levels and colors
levs = 0:0.005:.4;
cmap = colormap(flipud((brewermap(length(levs),'PuBuGn')))); %requires package 'brewermap'
%plot
contourf(KNramp,KSramp,AC_slope',levs,'EdgeAlpha',0);
hold on
cc = colorbar;
hold on
%add DH=0 boundary
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4]) 

xlim([2 50])
ylim([2 50])
clim([0.1 .4])
xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'$\frac{\delta\rho}{\delta\overline{\mathsf{H}}} (\mathsf{Sv}^{\mathsf{-1}})$','Interpreter','latex','FontSize',12)

%% normalized slopes
% for Figure S7 in the S.I. of Zimmerman et al (2024)

figure(11);clf
%set colormap levels and colors
levs = 0.3:0.02:1.5;
cmap = colormap(flipud((brewermap(length(levs)-1,'PuBuGn')))); %requires package 'brewermap'
%plot
contourf(KNramp,KSramp,var_norm',levs,'EdgeAlpha',0);
clim([0.4 1.5])
hold on
cc = colorbar;
hold on
%add DH=0 boundary
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4]) 

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
xlim([2 50])
ylim([2 50])
ylabel(cc,'$\frac{\delta\sigma^\mathsf{2}}{\delta\overline{\mathsf{H}}} (\mathsf{Sv}^{\mathsf{-1}})$','Interpreter','latex','FontSize',12)

figure(20);clf
%set colormap levels and colors
levs = 0.4:0.01:.7;
cmap = colormap(flipud((brewermap(length(levs),'PuBuGn')))); %requires package 'brewermap'
%plot
contourf(KNramp,KSramp,AC_norm',levs,'EdgeAlpha',0);
hold on
cc = colorbar;
hold on
%add DH=0 boundary
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4]) 

xlim([2 50])
ylim([2 50])
clim([0.4 .7])
xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'$\frac{\delta\rho}{\delta\overline{\mathsf{H}}} (\mathsf{Sv}^{\mathsf{-1}})$','Interpreter','latex','FontSize',12)


