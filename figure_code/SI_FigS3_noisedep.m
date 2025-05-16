% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Figure S3a and S3b in Zimmerman et al. 2024. Plots analytic solution and 
% stochastic integrations, 5 realizations and the average over all 
% realizations, of AMOC strength (Q) in the three box model for varying 
% hosing (H). Deterministic integrations conducted for various hosing 
% rates. Default calibration to FAMOUS_B 1xCO2.
%
% Dependencies:
%   - analytic_sol_3box_1xCO2_weak.mat or analytic_sol_3box_1xCO2_strong.mat 
%       (analytic equilibrium solutions of AMOC strength Q as a function of 
%       hosing H in the 3-box model for weak or strong gyres and model 
%       calibration FamousB_1xCO2)
%           -> different model calibrations (FamousB_2xCO2,HGEM,...) can be 
%           set here and will run in Analytic_3box.m 
%   - Qsol4_mathematica.mat (Unstable Q<0 analytic equilibrium solution of 
%       AMOC strength Q as a function of hosing H in the 3-box model for
%       calibration to FAMOUS_B 1xCO2. Solution derived from Mathematica)
%   - stochastic_3box_FMSB_1CO2_weak.mat OR
%       stochastic_3box_FMSB_1CO2_strong.mat (stochastic integrations of 
%       AMOC strength (Q) as a function of hosing (H) in the 3-box model 
%       for weak or strong gyres and model calibration FamousB_1xCO2) 
%           -> different model calibrations (FamousB_2xCO2,HGEM,...) can be 
%           set here and will run in Stochastic_sims_3box.m            
% 
% Output:
%   - Figure S3a and S3b in Zimmerman et al. 2024
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------
%% 
%choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
model = 'FMSB';
%choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
CO2 = 1;


% set gyre strength (1 = weak -> KN = 5.456 Sv, KS = 5.447 Sv; 2 = strong -> KN,KS = 27 Sv)
gyre_strength = 2;
gyre = {'weak', 'strong'};

%% load and plot analytic solution
% first load in or solve for the analytic solution associated with the 
% model calibration of choice 'analytic_sol_3box_FMSB_1CO2_weak.mat' and 
% 'analytic_sol_3box_FMSB_1CO2_strong.mat' available in Zenodo in 'Output' 
% folder
analytic_name = sprintf('analytic_sol_3box_%s_%dCO2_%s.mat',model,CO2,gyre{gyre_strength});
if exist(analytic_name,'file')==0
    Analytic_3box %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_name)

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength); clf
    hold on
    
    index_usn = find(h1>usn,1,'first'); %finds the index of the upper saddle node bifurcation
    
    plot(h1(1:index_usn),Qsol1(1:index_usn),'k','LineWidth',1.5) %plot positive stable solution from beginning until the usn
    hold on

    if gyre_strength == 1 %only weak gyre has unstable solution
        index_lsn = find(h1>lsn,1,'first'); %finds the index of the lower saddle node bifurcation
        index_zero = find(Qsol2<0,1,'last'); %finds the index where the unstable solution crosses Q=0

        plot(h1(index_zero:index_usn),Qsol2(index_zero:index_usn),'k--','LineWidth',1.5) %plot positive unstable solution from Q=0 until the upper saddle node
    else
        index_lsn = index_usn; %in the strong gyre scenario there is no bifurcation, this identifies the point where Q=0 and the solution switches from positive to negative
    end

    plot(h1(index_lsn:end),Qsol3(index_lsn:end),'k','LineWidth',1.5,'HandleVisibility','off') %plot negative stable solution

    ylabel('Q (Sv)','FontSize',14)
    xlabel('H(Sv)','FontSize', 15)
    
    xlim([-.2 0.3])
    ylim([-10 25])
    
    hold on
    box on

%%
%import data from mathematica for negative unstable solution (run for weak
%gyre only)
%Note: this solution is only valid for calibration to FAMOUS_B 1xCO2
load Qsol4_mathematica
H = squeeze(Expression1(1,:,:));
Qan = squeeze(Expression1(2,:,:));

figure(1);hold on;
plot(H,Qan,'k--','LineWidth',1.5,'HandleVisibility','off') %plot negative unstable solution

%% load and plot stochastic simulation runs
% load in or run the stochastic simulations associated with the 
% model calibration of choice 'stochastic_3box_FMSB_1CO2_weak.mat' and 
% 'stochastic_3box_FMSB_1CO2_strong.mat' available in Zenodo in 'Output' 
% folder
stoch_name = sprintf('stochastic_3box_%s_%dCO2_%s.mat',model,CO2,gyre{gyre_strength});
if exist(stoch_name,'file')==0
    Stochastic_sims_3box %If output does not exist, takes ~3 sec to run on a standard laptop Nov'24
end
load(stoch_name)

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength);
for l=1:2 %forward and backward
    Q_index = randsample(1000,5); %randomly select 5 realizations
    for i = 1:length(Q_index)
        plot(h_save(l,:),squeeze(Q_save(l,Q_index(i),:)),'linewidth',1,'Col',[.5 .5 1 .5])%plot 5 realizations
        hold on
    end
    plot(h_save(l,:),Qvec(l,:),'r','linewidth',1)%plot average over all realizations
end
%% Delta H*
%Only relevant for weak gyre scenario
%define Q value where hysteresis width is measured
thresh = 3;
%select hosing (Hbar) value associated with average crossing of above
%defined threshhold
husn = h_save(1,find(Qvec(1,:)<thresh,1,"first"));%forward
hlsn = h_save(2,find(Qvec(2,:)>thresh,1,"first"));%backward

%plot only on weak gyre plot
figure(1); hold on
plot([hlsn+.01 husn-.01],[thresh thresh],'-','linewidth',1,'Color',[0.4940 0.1840 0.5560])%line spanning hysteresis width (Delta H*)
plot([-0.2 0.3],[thresh thresh],'--','linewidth',.75,'Color',[0 0 0 .75]) %dashed line indicating threshhold level
%add arrows indicating span of DeltaH*
plot(hlsn+.008,thresh,'<','linewidth',1,'Color',[0.4940 0.1840 0.5560])
plot(husn-.008,thresh,'>','linewidth',1,'Color',[0.4940 0.1840 0.5560])
txt = '\Delta H*';
text(0,5,txt,'FontSize', 12)%add text label
box on