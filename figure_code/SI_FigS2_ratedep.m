% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Figure S2a and S2b in Zimmerman et al. 2024. Plots analytic solution and 
% deterministic integrations of AMOC strength (Q) in the three box model for
% varying hosing (H). Deterministic integrations conducted for various 
% hosing rates. Default calibration to FAMOUS_B 1xCO2.
%
% Dependencies:
%   - analytic_sol_3box_1xCO2_weak.mat or analytic_sol_3box_1xCO2_strong.mat 
%       (analytic equilibrium solutions of AMOC strength Q as a function of 
%       hosing H in the 3-box model for weak or strong gyres and model 
%       calibration FamousB_1xCO2)
%           -> different model calibrations (FamousB_2xCO2,HGEM,...) can be 
%           set here and will run in Analytic_3box.m 
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - solve_initial_salinities.m (solves for the steady state salinities 
%       (SN,ST) for initial hosing value (hmin) and final hosing value
%       (hmax) for the 3-box model)
%   - Qsol4_mathematica.mat (Unstable Q<0 analytic equilibrium solution of 
%       AMOC strength Q as a function of hosing H in the 3-box model for
%       calibration to FAMOUS_B 1xCO2. Solution derived from Mathematica)
%          
% Output:
%   - Figure S2a and S2b in Zimmerman et al. 2024
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
    
    xlim([-.5 0.5])
    
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

%%
%now run the 3-box model forward and backward with three different hosing
%rates and no noise

%load in parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
eval(params)

%set gyre strengths based on selection above
if gyre_strength == 1
    KN = 5.456; %Sv
    KS = 5.447;
    SNusn = 0.03446; %SN at the usn
    usn = 0.2138; % Hosing level when Q --> 0
else
    if gyre_strength == 2
        KN = 27;
        KS = 27;
        SNusn = 0.03446;
        usn = 0.2282;  % Hosing level when Q --> 0
    end
end

%set timestepping
dt = 1*syr; %x yr in seconds
%set hosing rates
rate = [.0005,.00005,.00001]; %Sv/year
rate_legend = {sprintf('%0.0e Sv/yr',rate(1)),sprintf('%g Sv/yr',rate(2)),sprintf('%g Sv/yr',rate(3))};

cc = ['g','m','c']; %plot colors for each rate

%hosing range
hmax = .5;
hmin = -1;

%get initial salinities (SN,ST) for hmin and hmax
solve_initial_salinities

for k = 1:length(rate) %run each hosing rate
    %set rate 
    j = rate(k); %seconds
    tend = (hmax-hmin)/j; %years
    runtime = tend*syr; %seconds
    nt = round(runtime/dt);
    t = linspace(0,tend,nt); %years
    h1 = linspace(hmin,hmax,nt); %hosing vector forward
    h2 = linspace(hmax,hmin,nt); %hosing vector backward

        for l = 1:2 %forward and backward
            %load in initial conditions
            SN = SNsol0(l);
            ST = STsol0(l);
            
            Qvec = nan(nt,1); %resets each run, outputs not saved here, only plotted
            %select hosing for forward or backward run
            if l == 1
                h = h1;
            else
                h = h2;
            end
            tic
            for i = 1:nt %timestepping
                H = h(i);
                SIP = 1/VIP*(C-(VN*SN+VT*ST+VS*SS+VD*SD)); %salt conservation
                Q = ls*(a*(TS-T0)+b*(SN-SS));
        
                if Q>0
                    dSN = 1/VN*(Q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0);
                    dST = 1/VT*(Q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0);
                elseif Q<0
                    dSN = 1/VN*(abs(Q)*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0);
                    dST = 1/VT*(abs(Q)*(SN-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0);
                end
                %update salinities
                SN = SN + dSN*dt;
                ST = ST + dST*dt;
                Qvec(i) = Q; %save for plot
            end
            toc
            %figure(1) will be weak gyre scenario; figure(2) will be strong 
            % gyre scenario
            figure(gyre_strength) 
            if l==1
                plot(h,Qvec,'-','col',cc(k),'LineWidth',.75) %plot forward run
            else
                plot(h,Qvec,'-','col',cc(k),'LineWidth',.75,'HandleVisibility','off') %plot backward run; want each rate in legend only once
            end
            hold on
            axis tight 
            xlim([-0.5 0.5])
        end
end
figure(1); %legend for weak gyres only; if you want legend for strong gyre, remove 'Unstable solution'
legend('Stable solution','Unstable solution',rate_legend{1},rate_legend{2},rate_legend{3})