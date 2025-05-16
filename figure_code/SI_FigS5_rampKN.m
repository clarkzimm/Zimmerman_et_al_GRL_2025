% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Figure S5 in Zimmerman et al. 2024. Plots analytic solution and 
% stochastic integrations of AMOC strength (Q) in the three box model for
% varying North Atlantic gyre transport strength (KN). calibrated to 
% FAMOUS_B 1xCO2. Computes and plots statistical indicators variance, 
% autocorrelation, and de-correlated autocorrelation (Boers 2021) for 
% plotted stochastic runs.
%
% Dependencies:
%   - Analytic_3box_rampKN.m (solves for the analytic equilibrium solutions 
%       of AMOC strength Q as a function of North Atlantic gyre transport 
%       strength (KN) in the 3-box model for various model calibrations 
%       (FamousB_1xCO2,HGEM,...)
%   - Stochastic_sims_3box_rampKN.m (runs stochastic integrations of 
%       AMOC strength (Q) as a function of North Atlantic gyre transport 
%       strength (KN) in the 3-box model for various model calibrations 
%       (FamousB_1xCO2,HGEM,...)   
%          
% Output:
%   - Figure S5 in Zimmerman et al. 2024
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

%% load and plot analytic solution
% first load in or solve for the analytic solution associated with the 
% model calibration of choice 'analytic_sol_3box_FMSB_1CO2_rampKN.mat' 
% available in Zenodo in 'Output' folder
analytic_name = sprintf('analytic_sol_3box_%s_%dCO2_rampKN.mat',model,CO2);
if exist(analytic_name,'file')==0
    Analytic_3box_rampKN %If output does not exist, takes ~3 min to run on a standard laptop Nov'24
end
load(analytic_name)

figure(1); clf
    hold on
    subplot(2,2,[1 3]); hold on %selects first column for plot
    
    plot(kn1,Qsol1,'k','LineWidth',1.5) %plot positive stable solution
    hold on

    ylabel('Q (Sv)','FontSize',14)
    xlabel('$K_N (Sv)$','FontSize',15,'Interpreter','latex')
    set ( gca, 'xdir', 'reverse' ) %plot with decreasing gyre strength to the right
    ylim([10 20])
    
    hold on

%% load and plot stochastic simulations
% load in or run the stochastic simulations associated with the 
% model calibration of choice 'stochastic_3box_FMSB_1CO2_rampKN.mat' 
% available in Zenodo in 'Output' folder
stoch_name = sprintf('stochastic_3box_%s_%dCO2_rampKN.mat',model,CO2);
if exist(stoch_name,'file')==0
    Stochastic_sims_3box_rampKN %If output does not exist, takes ~1 sec to run on a standard laptop Nov'24
end
load(stoch_name)

figure(1);
subplot(2,2,[1 3]); hold on %plots on top of the analytic solution
Q_index = randsample(1000,5); %randomly selects 5 realizations to plot
for i = 1:length(Q_index)
    plot(KN_save(2,:),squeeze(Q_save(2,Q_index(i),:)),'linewidth',1,'Col',[.5 .5 1 .5]) %plots backward runs only (want decreasing KN)
    hold on
end
ax = gca;
set(ax,'Ticklength',[0.01 0.025], 'FontSize', 14)
set ( gca, 'xdir', 'reverse' )

%% time axis
axis tight
box on
hold on

figure(1);hold on; subplot(2,2,[1 3]); hold on
hAx =gca;
hAx(2) = axes('position',hAx(1).Position, 'color','none','XAxisLocation','top','YAxisLocation','right','FontSize',14); %selects top x axis
xlim([0 time(end)])
xticks(0:1000:time(end))  %sets ticks for time axis
xlabel('Time (yrs)','FontSize', 14)
yticks([])
set(hAx, 'Position', [.05, 0.15, .5, .75]);

%% statistics
% now statistical indicators including decorrelated AC
% first get anomaly
Qvec_re = squeeze(Q_save(2,:,:)); %select backward run
Qdet = squeeze(Qdet(2,:)); %select backward run
int = 10; %define lag
q_anom = Qvec_re-Qdet; %get anomaly
q_anom_short = q_anom(:,1:int:end); %cut vector for lagged autocorrelation

pick = Q_index; %select same five realizations as plotted above for calculating statistics

% time window
w = round(.05*length(q_anom_short(1,:))); %note this is half the window
kn3 = kn1(1:int:L+int); %hosing with lag accounted for
KN_ind = kn3(w:end-w); %hosing indexed at center of running window

varian = nan(1,length(KN_ind)); rho = varian; rhob = varian;
varian2 = varian; rho2 = varian2; rhob2 = varian2;
tic
for i = 1+w:length(KN_ind)+w-3
    q1 = q_anom_short(:,i+(-w:w)); %selects each window for all realizations
    q2 = q_anom_short(:,i+(-w:w)+1); %window one lag behind q1
    q3 = q_anom_short(:,i+(-w:w)+2); %window one lag behind q2
    q4 = q1(pick,:); %selects each window for 5 realizations
    q5 = q2(pick,:);  %window one lag behind q4
    q6 = q3(pick,:); %window one lag behind q5
    %first variance
    varian(i-w) = var(q2(:)); %variance of all realizations
    varian2(i-w) = var(q5(:)); %variance of 5 selected realizations
    %then normal correlation
    p = corrcoef(q2(:),q1(:)); %correlation coefficients of all realizations
    p2 = corrcoef(q5(:),q4(:)); %correlation coefficients of 5 realizations
    rho(i-w) = p(1,2); %selects ac
    rho2(i-w) = p2(1,2);%selects ac
    %then decorrelated autocorrelation
    V1 = q3 - p(1,2).*q2; %follows from Boettner and Boers (2022) equation 16
    V2 = q2 - p(1,2).*q1;
    V3 = q6 - p2(1,2).*q5;
    V4 = q5 - p2(1,2).*q4;
    pb = corrcoef(V1(:),V2(:));
    pb2 = corrcoef(V3(:),V4(:));
    rhob(i-w) = pb(1,2);
    rhob2(i-w) = pb2(1,2);
end
toc

%find unbiased estimator from decorrelated autocorrelation rhob
%follows from Boettner and Boers (2022) equation 18
%for all realizations
A = rho+rhob; B = rhob./rho;
rhoa = 1/2*(A+sqrt(A.^2-4*B));
%and for 5 realizations
A2 = rho2+rhob2; B2 = rhob2./rho2;
rho2a = 1/2*(A2+sqrt(A2.^2-4*B2));

%% plot

figure(1); 

% plot variance
aa = subplot(2,2,2);%top right subplot
plot(KN_ind,varian,'b','LineWidth',1.5) %plot full variance
hold on
plot(KN_ind,varian2,'LineWidth',1,'Col',[.5 .5 1 .75]) %plot 5 realization variance
ylim([0 1]) 

xticklabels([])
box on
set(aa, 'Position', [.62, 0.55, .31, .35]);
set(gca,'Ticklength',[0.01 0.025], 'FontSize', 14)
set (gca, 'xdir', 'reverse' )
ylabel('\sigma^2','FontSize',16)


% plot autocorrelation
bb = subplot(2,2,4); %bottom right subplot
plot(KN_ind,rho,'b','LineWidth',1.5) %plot full AC
hold on
plot(KN_ind,rho2,'LineWidth',1,'Col',[.5 .5 1 .75],'LineStyle','-') %plot 5 realizations ac
ylim([rho(plot_min)-.005 rho(end-3)+.005]) %y limit just above and below data range
plot(KN_ind,rhoa,'Col',[0 0.7 0],'LineWidth',1.5) %plot full unbiased estimator
plot(KN_ind,rho2a,'LineWidth',1,'Col',[0 1 0 .5],'LineStyle','-') %plot 5 realizations estimator

ylim([0.4 1])
set(bb, 'Position', [.62, 0.15, .35, .35]);
set(gca,'Ticklength',[0.01 0.025], 'FontSize', 14)
set (gca, 'xdir', 'reverse' )  %plot with decreasing gyre strength to the right

ylabel('\rho,\rho_a','FontSize',16)
xlabel('$K_N (Sv)$','FontSize',15,'Interpreter','latex')




