% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Creates Figure 4 from Zimmerman et al (2024). Plots the width of the
% bistable solution (DH) for the 5-box model, the difference between DH3 
% and DH5, the width of the bistable solution for the 3-box GW calibration, 
% and the difference between DH1xCO2 and DH2xCO2 as functions of transport 
% strengths (KN,KS). For Figure 4, all for model calibrations to FAMOUS_B. 
% To create Figure S8, use model calibration for HadGEM-AO.
%
% Inputs:
%   - analytic_DH_Xbox_YY_ZCO2.mat (analytic equilibrium hysteresis (DH) in 
%       hosing (H) for varying gyre strengths (KN,KS), and boundary between 
%       monostable and bistable (KN,KS) regimes for X box model, with YY 
%       ZxCO2 parameterization)
%
% Dependancies:
%   - Analytic_DeltaH.m (solves for analytic equilibrium hysteresis (DH) in 
%       hosing (H) for varying gyre strengths (KN,KS), and boundary between 
%       monostable and bistable (KN,KS) regimes if output file for
%       model/param choice does not exist yet)
%          
% Output:
%   - Figure 4a,b,c and d from Zimmerman et al (2024)
%   - Figure S8 from Zimmerman et al (2024)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
%close all; clear all;

%% Choose Model run
% for Figure 4 in Zimmerman et al (2024) set model = 'FMSB'
% for Figure S8 in the S.I. of Zimmerman et al (2024) set model = 'HGEM'
%choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
model = 'FMSB';

%choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW)) Note: only
%adjusts 3- vs 5- box comparison (1 in Zimmerman et al (2024))
CO2 = 1;
%% 5- vs 3- box
% first load in or solve for the width of the unstable analytic solution 
% associated with the model calibration of choice for 5 box formulation
% 'analytic_DH_5box_FMSB_1CO2.mat' available in Zenodo in 'Output' folder
analytic_3 = sprintf('analytic_DH_3box_%s_%dCO2.mat',model,CO2);
analytic_5 = sprintf('analytic_DH_5box_%s_%dCO2.mat',model,CO2);
if exist(analytic_5,'file')==0
    box = 5; %set to 5-box formulation
    Analytic_DeltaH %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_5)
%rename for comparison
DH5 = DH;
DHbound5 = DHbound;
%load in 3-box DHbound 'analytic_DH_3box_FMSB_1CO2.mat' available in Zenodo 
% in 'Output' folder
if exist(analytic_3,'file')==0
    box = 3; %set to 3 box formulation
    Analytic_DeltaH %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_3)

%First plot the 5-box model (Figure 4a in Zimmerman et al (2024))
figure(28); clf
%define colormap levels and colors
levs = 0:0.025:.4;
cmap = colormap((brewermap(length(levs)-1,'YlGnBu'))); %NOTE: requires package 'brewermap'
contourf(KNvec,KSvec,DH5',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
hold on

cc = colorbar;
% add DH=0 bounds for 3- and 5-box formulations
plot(KNvec,DHbound,'k--','linewidth',1)
plot(KNvec,DHbound5,'k-','linewidth',2)

clim([0 0.4])
xlim([0 50])
ylim([0 50])

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH (Sv)')
%%
%add markers for different calibrations if desired
hold on
%Gyre strength values for various EMICs from Wood et al. 2015
KN_GCM = [5.439 1.762 5.601 15.890 20.954];
KS_GCM = [1.880  1.872 7.169 6.828 8.384];
%define labels
GCMs = {'3-box bound','5-box bound','Scenario 1','Scenario 2','FMS_A 1xCO_2', 'FMS_B 2xCO_2' ,'HG 1xCO_2', 'HG 2xCO_2' ,'HG 4xCO_2'};
%define marker colors
c = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560], [0 0.4470 0.7410]};
%markers for weak and strong gyre scenarios
scatter(5.456,5.447,100,'m','hexagram','filled')
scatter(27,27,100,'c','hexagram','filled')
%plot markers
for i=1:5
    scatter(KN_GCM(i),KS_GCM(i),50,cell2mat(c(i)),"filled")
end
%add legend
legend(GCMs)
hold  off
%%
%now get difference between 3- and 5-box formulations
DDH = DH-DH5;
zeromask = DDH < 1e-8 & DDH > -1e-8; 
DDH(zeromask) = 0; %set numerical differences to zero

%set colormap
levs = -.5:0.0001:.5;
col = colormap((brewermap(length(levs),'PRGn'))); %NOTE: requires package 'brewermap'

%plot Figure 4b in Zimmerman et al (2024)
figure(34); clf
contourf(KNvec,KSvec,DDH',levs,'EdgeAlpha',0,'HandleVisibility','off')
colormap(col)
cc = colorbar;
clim([-.04 .04])
hold on
%add DH=0 bounds
plot(KNvec,DHbound,'k--','linewidth',1)
plot(KNvec,DHbound5,'k-','linewidth',1)
legend('3-box bound','5-box bound') 

xlim([0 50])
ylim([0 50])
xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\Delta H(3-box) - \Delta H(5-box) (Sv)')
%%
%add additional labels
txt1 = '\leftarrow';% 3-box more bistable';
txt2 = '\leftarrow';% 5-box more bistable';
text(8,5,txt2,'Color',[0.4940 0.1840 0.5560],'FontSize',20)
text(19,36,txt1,'Color',[0.4660 0.6740 0.1880],'FontSize',20)
%% 1xCO2 vs 2xCO2
% set model formulation(3 -> 3-box; 5 -> 5-box) (3 in Zimmerman et al (2024))
box = 3;

% load in or solve for the width of the unstable analytic solution 
% associated with the model formulation of choice for 2xCO2 calibration
% 'analytic_DH_3box_FMSB_2CO2.mat' and 'analytic_DH_3box_HGEM_2CO2.mat' 
% available in Zenodo in 'Output' folder
analytic_2x = sprintf('analytic_DH_%dbox_%s_2CO2.mat',box,model);
analytic_1x = sprintf('analytic_DH_%dbox_%s_1CO2.mat',box,model);
if exist(analytic_2x,'file')==0
    CO2 = 2;%set 2xCO2
    Analytic_DeltaH %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_2x)
%rename for comparison
DH2x = DH;
DHbound2x = DHbound;
% load in 1xCO2 DHbound 'analytic_DH_3box_FMSB_1CO2.mat' and 
% 'analytic_DH_3box_HGEM_1CO2.mat' available in Zenodo in 'Output' folder
if exist(analytic_1x,'file')==0
    CO2 = 1;%set 1xCO2
    Analytic_DeltaH
end
load(analytic_1x)

%plot Figure 4c in Zimmerman et al (2024) if model = 'FMSB'
%plot Figure S8b in Zimmerman et al (2024) if model = 'HGEM'
figure(21);clf
%define colormap
levs = 0:0.025:.4;
cmap = colormap((brewermap(length(levs)-1,'YlGnBu')));%NOTE: requires package 'brewermap'
contourf(KNvec,KSvec,DH2x',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
hold on
cc = colorbar;
%add DH=0 bounds
plot(KNvec,DHbound2x,'k:','linewidth',2,'HandleVisibility','off')
plot(KNvec,DHbound,'k--','linewidth',1,'HandleVisibility','off')
clim([0 0.4])

xlim([0 50])
ylim([0 50])

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH (Sv)')
%%
%now get difference
DDH = DH-DH2x;
zeromask = DDH < 1e-8 & DDH > -1e-8;
DDH(zeromask) = 0; %set numerical differences to zero
%choose colormap
levs = 0:0.0001:.2;
col = colormap((brewermap(length(levs),'Greens')));%NOTE: requires package 'brewermap'
col(1,:) = [1 1 1]; %white out diff = 0

%plot Figure 4d in Zimmerman et al (2024)
figure(33); clf
contourf(KNvec,KSvec,DDH',levs,'EdgeAlpha',0,'HandleVisibility','off')
colormap(col)
cc = colorbar;
hold on
%add DH=0 bounds
plot(KNvec,DHbound,'k--','linewidth',1)
plot(KNvec,DHbound2x,'k:','linewidth',1)
legend('1xCO_2 bound','2xCO_2 bound') 

xlim([0 50])
ylim([0 50])
xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\Delta H(1xCO_2) - \Delta H(2xCO_2) (Sv)')
%%
%plot Figure S8a in the S.I. of Zimmerman et al (2024)
figure(22);clf
%define color map
levs = 0:0.025:.4;
cmap = colormap((brewermap(length(levs)-1,'YlGnBu')));%NOTE: requires package 'brewermap'
%plot
contourf(KNvec,KSvec,DH',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
hold on
cc = colorbar;
%add DH=0 bound
plot(KNvec,DHbound,'k--','linewidth',2,'HandleVisibility','off')
clim([0 0.4])

xlim([0 50])
ylim([0 50])

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH (Sv)')