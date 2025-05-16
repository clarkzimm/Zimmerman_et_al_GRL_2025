% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Creates Figure S6a and S6b from Zimmerman et al (2024). Plots the width 
% of the hysteresis (DH*) experienced by deterministic simulations of the 
% 3-box model, and the difference between DH* and the analytic width of the
% bistable solution (DH). For Figure S6, model calibrated to FAMOUS_B 1xCO2.
%
% Dependancies:
%   - analytic_DH_3box_FMSB_1CO2_short.mat (analytic equilibrium hysteresis (DH) 
%       in hosing (H) for varying gyre strengths (KN,KS), and boundary 
%       between monostable and bistable (KN,KS) regimes for 3 box model, 
%       calibrated to FAMOUS_B 1xCO2)
%           - Run 'Analytic_DeltaH.m' if alternate model calibration is
%                   desired. Set 'res = 100' there so density matches
%                   deterministic data.
%   - Det_D_3box.m (solves for the width of the hysteresis in hosing (H) 
%       experienced by deterministic runs of the 3-box model for varying 
%       gyre strengths (KN,KS))
%          
% Output:
%   - Figure S6a and S6b from Zimmerman et al (2024)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
%% S6 a
%choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
model = 'FMSB';
%choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
CO2 = 1;

% load in or solve for the deterministic hysteresis width for the model 
% calibration of choice 'det_DH_3box_FMSB_1CO2.mat' available in Zenodo in 
% 'Output' folder
det_name = sprintf('det_DH_3box_%s_%dCO2.mat',model,CO2);
if exist(det_name,'file')==0
    Det_DH_3box %If output does not exist, takes ~25 sec to run on a standard laptop Nov'24
end
load(det_name)

%load in the analytic width of the bistable solution, available in Zenodo 
% in 'Output' folder
load analytic_DH_3box_FMSB_1CO2_short.mat %if alternate parameterization is desired, run Analytic_DeltaH with model calibration of choice and res=100

figure(23);clf
%define colormap levels and colors
levs = 0:0.033:.4;
cmap = colormap((brewermap(length(levs),'YlGnBu'))); %NOTE: requires package 'brewermap'
%plot figure S6a in Zimmerman et al (20024)
contourf(KNramp,KSramp,detDH',levs,'EdgeAlpha',0.25,'HandleVisibility','off')
hold on
cc = colorbar; set(cc,'ticklabels',levs(1:2:end),'ticks',levs(1:2:end))
%add DH=0 bound
plot(KNvec,DHbound,'--','Linewidth',2,'Color',[0 0 0 .4])
clim([0 0.4])

xlim([2 50])
ylim([2 50])

xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\DeltaH^* (Sv)')

%% S6 b
%compute difference between analytic width of bistable solution and
%deterministic hysteresis width
DDH = detDH-DH;

%set colormap levels and colors
levs = 0:0.001:.2;
col = colormap((brewermap(length(levs),'Greens'))); %NOTE: requires package 'brewermap'
col(1,:) = [1 1 1]; %set 0 difference to white

%plot Figure S6b in Zimmerman et al (2024)
figure(32); clf
contourf(KNvec,KSvec,DDH',levs,'EdgeAlpha',0,'HandleVisibility','off')
colormap(col)
cc = colorbar;
hold on
%add DH=0 boundary
plot(KNvec,DHbound,'k--','linewidth',1)

clim([0 0.1])
xlim([2 50])
ylim([2 50])
xlabel('K_N (Sv)')
ylabel('K_S (Sv)')
ylabel(cc,'\Delta H^* - \Delta H (Sv)')