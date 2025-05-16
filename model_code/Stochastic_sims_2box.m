% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Integrates the 2-box model for both forward and backward hosing (H),
% solving for AMOC strength (Q), with Guassian white noise added to the
% hosing. Called by SI_FigS4_2box.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - solve_initial_salinities_2box.m (solves for the steady state 
%       salinities (SN,ST) for initial hosing value (hmin) and final hosing 
%       value (hmax)
%          
% Output:
%   -'h_save','h1','h2' (deterministic hosing (bar(H)), h1- forward, h2-
%       backward, h_save = [h1 h2]
%   -'Qvec(2,nt)' (noisy simulations meaned over the realizations, 1 - 
%       forward, 2 - backward) 
%   -'Q_save(2,re,nt)' (full noisy simulations for all realizations)
%   -'Qdet(2,nt)' (deterministic simulations)
%   -'H_bif_min' (index of first crossing of Qusn, defined below (default =
%       0)
%   -'time' (time vector for ploting)    
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
%close all; clear all;

tic
% %uncomment if running here, comment if calling from other script
% %NOTE: usn, lsn, and hopf hosing values are only valid for FAMOUS_B 1xco2
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM2-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;
%
% % set gyre strength (1 = weak -> KN = 5.456 Sv, KS = 5.447 Sv; 2 = strong -> KN,KS = 27 Sv)
% gyre_strength = 2;

%load in parameters for callibration defined above or in SI_FigS4_2box.m
params = sprintf('parameters_3box_%s_%dCO2',model,CO2); %note: same parameter file as 3-box; almost all parameters stay the same
eval(params)
SIP = SIP0; %Now SIP is held constant

if gyre_strength == 1
    KN = 5.456; %Sv
    %set hosing range; different for weak/strong
    hmax = 1.5; % ending hosing values
    hmin = -2; %initial hosing value
else
    if gyre_strength == 2
        KN = 22;%Sv
        %set hosing range; different for weak/strong
        hmax = 10.5; % ending hosing values
        hmin = 7; %initial hosing value
    end
end

%define timestep and hosing
dt = 1*syr; %x yr in seconds
rate = .0005; %Hosing rate in Sv/year

tend = (hmax-hmin)/rate; %runtime in years
runtime = tend*syr; %seconds
nt = round(runtime/dt);
time = linspace(0,6*tend/7,6*nt/7); %for plot (only plotting 6/7th of total integration period) in years

h1 = linspace(hmin,hmax,nt); %forward hosing deterministic
h2 = linspace(hmax,hmin,nt); %backward hosing deterministic

%add noise
sig1 = 3e4; %set noise level
sig = sig1/sqrt(dt); %normalize by timestep
re = 1000; %set # of realizations

N1 = sig*randn(re,nt); %Gaussian white noise with var sig^2
N2 = sig*randn(re,nt);
H1 = repmat(h1,re,1)+N1; %forward hosing noisy
H2 = repmat(h2,re,1)+N2; %backward hosing noisy

Qusn = 0; % define shut down

%get initial conditions
solve_initial_salinities_2box

clear h_save
%define/save variables
h_save(1,:) = h1;
h_save(2,:) = h2;
Q_save = nan(2,re,nt);
Qvec = nan(2,nt);
Qdet = nan(2,nt);
H_bif = nan(re,1);

for k = 1:2 %forward and backward
    Qvec_re = nan(re,nt);
    %set hosing
    if k == 1 %forward
        h = h1; %Hbar
        HH = H1; %H
    else %backward 
        h = h2; %Hbar
        HH = H2; %H
    end
    for p =1:re+1 %step through realizations and last run is deterministic
        SN = SNsol0(k); %load initial salinity, resets every realization

        for i = 1:nt %timestepping
            if p==re+1 %last run uses deterministic hosing Hbar
                H = h(i); %deterministic hosing
            else
                H = HH(p,i); %noisy hosing
            end
            
            ST = (C-(VN*SN+VIP*SIP+VS*SS+VD*SD))/VT; %salt conservation
            Q = ls*(a*(TS-T0)+b*(SN-SS));
    
            if Q>0
                dSN = 1/VN*(Q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0);
            elseif Q<0
                dSN = 1/VN*(abs(Q)*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0);
            end
            %update salinity
            SN = SN + dSN*dt;
    
            if p == re+1
                Qdet(k,i) = Q; %save det run
            else
                Qvec_re(p,i) = Q; %save stochastic runs
            end
        end 
        
        if p < re+1 %only for stochastic runs with increasing hosing
            if k == 1
                H_bif(p) = find(Qvec_re(p,:)<Qusn,1,"first"); %index of shut down
            end
        end
    end
    
    if gyre_strength == 1
        H_bif_min = min(H_bif); %in the weak gyre case: index of earliest bifurcation
    else
        if gyre_strength == 2
            H_bif_min = mean(H_bif); %in the strong gyre case: no bifurcation, so, mean index of shut down
        end
    end
    
    Qvec(k,:) = mean(Qvec_re); %avg stochastic run
    Q_save(k,:,:) = Qvec_re; %save all runs
    
end
toc
%%
% save output
save(stoch_name,'h1','h2','h_save','Qvec','Q_save','Qdet','time','H_bif_min')
