% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Figure S4a and S4b in Zimmerman et al. 2024. Plots analytic solution and 
% stochastic integrations of AMOC strength (Q) in the two box model for
% varying hosing (H). calibrated to FAMOUS_B 1xCO2. Computes and plots 
% statistical indicators variance, autocorrelation, and de-correlated 
% autocorrelation (Boers 2021) for plotted stochastic runs.
%
% Dependencies:
%   - analytic_sol_2box_1xCO2_weak.mat or analytic_sol_2box_1xCO2_strong.mat 
%       (analytic equilibrium solutions of AMOC strength Q as a function of 
%       hosing H in the 2-box model for weak or strong gyres and model 
%       calibration FamousB_1xCO2)
%           -> different model calibrations (FamousB_2xCO2,HGEM,...) can be 
%           set here and will run in Analytic_2box.m 
%   - stochastic_2box_FMSB_1CO2_weak.mat OR
%       stochastic_2box_FMSB_1CO2_strong.mat (stochastic integrations of 
%       AMOC strength (Q) as a function of hosing (H) in the 2-box model 
%       for weak or strong gyres and model calibration FamousB_1xCO2) 
%           -> different model calibrations (FamousB_2xCO2,HGEM,...) can be 
%           set here and will run in Stochastic_sims_2box.m    
%          
% Output:
%   - Figure S4a and S4b in Zimmerman et al. 2024
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
% model calibration of choice 'analytic_sol_2box_FMSB_1CO2_weak.mat' and 
% 'analytic_sol_2box_FMSB_1CO2_strong.mat' available in Zenodo in 'Output' 
% folder
analytic_name = sprintf('analytic_sol_2box_%s_%dCO2_%s.mat',model,CO2,gyre{gyre_strength});
if exist(analytic_name,'file')==0
    Analytic_2box %If output does not exist, takes ~2 min to run on a standard laptop Nov'24
end
load(analytic_name)

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength); clf
    hold on
    subplot(2,2,[1 3]); hold on %selects first column for plot
    
    index_usn = find(han>usn,1,'first'); %finds the index of the upper saddle node bifurcation
    
    plot(han(1:index_usn),Qsol1(1:index_usn),'k','LineWidth',1.5) %plot positive stable solution from beginning until the usn
    hold on

    if gyre_strength == 1 %only weak gyre has unstable solution
        index_zero = find(han>lsn,1,'first'); %finds the index of the lower saddle node bifurcation; in 2-box model this is also the index where the negative solution crosses Q=0

        plot(han(index_zero:index_usn),Qsol2(index_zero:index_usn),'k--','LineWidth',1.5) %plot positive unstable solution from Q=0 until the upper saddle node
        Qusn = Qsol1(index_usn); %Finds Q at the upper saddle node
        Qlsn = Qsol3(index_zero); %Finds Q at the lower saddle node

        % stars at the saddle-node bifurcations 
        plot(usn,Qusn,'pentagram','MarkerSize',14,'MarkerFaceColor','g','MarkerEdgeColor','k')
        plot(lsn,Qlsn,'pentagram','MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor','k')
    else
        index_zero = index_usn; %in the strong gyre scenario there is no bifurcation, this identifies the point where Q=0 and the solution switches from positive to negative
    end

    plot(han(index_zero:end),Qsol3(index_zero:end),'k','LineWidth',1.5) %plot negative stable solution
    x1 = xlim;
    ylabel('Q (Sv)','FontSize',14)
    xlabel('$\overline{\mathsf{H}}$ $\mathsf{(Sv)}$','interpreter','latex', 'FontSize', 15)
    
    hold on
    box on

%% load and plot stochastic simulations
% first load in or run the stochastic simulations associated with the 
% model calibration of choice 'stochastic_2box_FMSB_1CO2_weak.mat' and 
% 'stochastic_2box_FMSB_1CO2_strong.mat' available in Zenodo in 'Output' 
% folder
stoch_name = sprintf('stochastic_2box_%s_%dCO2_%s.mat',model,CO2,gyre{gyre_strength});
if exist(stoch_name,'file')==0
    Stochastic_sims_2box %If output does not exist, takes ~1 sec to run on a standard laptop Nov'24
end
load(stoch_name)

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength);
subplot(2,2,[1 3]); hold on %plots on top of the analytic solution
Q_index = randsample(1000,5); %randomly selects 5 realizations to plot
for i = 1:length(Q_index)
    plot(h_save(1,:),squeeze(Q_save(1,Q_index(i),:)),'linewidth',1,'Col',[.5 .5 1 .5]) %plots forward runs only
    hold on
end
ax = gca;
ax.FontSize=14;
xlim(x1)
%% delta H label
% add delta H label (run for weak gyre only)
figure(1);hold on; subplot(2,2,[1 3]); hold on %plots on top of analytic and stochastic plot

plot([lsn+.01 usn-.01],[3 3],'k-','linewidth',1) %line from the lsn to the usn at Q=3
plot(lsn+.005,3,'k<','linewidth',1) %arrow pointing towards lsn
plot(usn-.005,3,'k>','linewidth',1) %arrow pointing towards usn
txt = '\Delta H';
text(-1,5,txt,'FontSize', 16) %adds text label
%% line of earliest shut off and time axis
%add vertical line at the hpsing (H) where the first realization goes to
%zero (Q-->0)
y1 = ylim; plot(h1(H_bif_min)*[1 1],y1,':','Color',[1 0 1],'linewidth',1)
axis tight 
box on
hold on

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength);hold on; subplot(2,2,[1 3]); hold on
hAx =gca;
hAx(2) = axes('position',hAx(1).Position, 'color','none','XAxisLocation','top','YAxisLocation','right','FontSize',14); %selects top x axis
xlim([0 time(end)]) 
xticks(0:1000:time(end)) %sets ticks for time axis
xlabel('Time (yrs)','FontSize', 14)
yticks([])
set(hAx, 'Position', [.05, 0.15, .5, .75]);

%% statistics
% now statistical indicators including decorrelated AC
% first get anomaly
Qvec_re = squeeze(Q_save(1,:,:)); %select forward run
Qdet = squeeze(Qdet(1,:)); %select forward run
L = round(H_bif_min); %index of first realization to reach Q=0
int = 10; %define lag
q_anom = Qvec_re-Qdet; %get anomaly
q_anom_short = q_anom(:,1:int:L+int); %cut vector for lagged autocorrelation and end at index of first shut off

pick = Q_index; %select same five realizations as plotted above for calculating statistics

% time window
w = round(.05*length(q_anom_short(1,:))); %note this is half the window
h3 = h1(1:int:L+int); %hosing with lag accounted for
H_ind = h3(w:end-w); %hosing indexed at center of running window

varian = nan(1,length(H_ind)); rho = varian; rhob = varian;
varian2 = varian; rho2 = varian2; rhob2 = varian2;
tic
for i = 1+w:length(H_ind)+w-3
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

%% plot statistics
bif_min = h1(L); %select hosing level of first crossing of Q=0
plot_min = find(H_ind>x1(1),1,'first'); %find index of lower x limit for plotting

%figure(1) will be weak gyre scenario; figure(2) will be strong gyre
%scenario
figure(gyre_strength); 

% plot variance
aa = subplot(2,2,2);%top right subplot
plot(H_ind,varian,'b','LineWidth',1.5) %plot full variance
hold on
plot(H_ind,varian2,'LineWidth',1,'Col',[.5 .5 1 .75]) %plot 5 realization variance
ylim([varian(plot_min)-.1 varian(end-3)+.2]) %y limit just above and below data range

y1 = ylim; plot(bif_min*[1 1],y1,':','Color',[1 0 1],'linewidth',1) %plot vertical line where first realization crosses Q=0
xticklabels([])
box on
set(aa, 'Position', [.62, 0.55, .31, .35]);
set(gca,'Ticklength',[0.01 0.025], 'FontSize', 14)
ylabel('\sigma^2','FontSize',16)

% plot autocorrelation
bb = subplot(2,2,4); %bottom right subplot
plot(H_ind,rho,'b','LineWidth',1.5) %plot full AC
hold on
plot(H_ind,rho2,'LineWidth',1,'Col',[.5 .5 1 .75],'LineStyle','-') %plot 5 realizations ac
plot(H_ind,rhoa,'Col',[0 0.7 0],'LineWidth',1.5) %plot full unbiased estimator
plot(H_ind,rho2a,'LineWidth',1,'Col',[0 .7 0 .5],'LineStyle','-') %plot 5 realizations estimator

ylim([rho(plot_min)-.005 rho(end-3)+.005]) %y limit just above and below data range
ylabel('\rho , \rho_{corr}','FontSize',16)
xlabel('$\overline{\mathsf{H}}$ $\mathsf{(Sv)}$','FontSize',15,'Interpreter','latex')

y1 = ylim; plot(bif_min*[1 1],y1,':','Color',[1 0 1],'linewidth',1) %plot vertical line where first realization crosses Q=0

% %uncomment to plot window size
% plot([bif_min-0.225 bif_min],[.92 .92],'k-|') 
% txt = 'W';
% text(bif_min-0.15,.935,txt)



