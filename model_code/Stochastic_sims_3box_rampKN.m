% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Integrates the 3-box model for both forward and backward ramping of North 
% Atlantic gyre transport strength (KN), solving for AMOC strength (Q), 
% with Guassian white noise added to KN. Called by SI_FigS5_rampKN.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - solve_initial_salinities_KN.m (solves for the steady state salinities 
%       (SN,ST) for initial gyre value (knmax) and final gyre value (knmin)
%          
% Output:
%   -'KN_save','kn1','kn2' (deterministic gyre ramping (bar(KN)), kn1- 
%       forward, kn2-backward, kn_save = [kn1 kn2]
%   -'Qvec(2,nt)' (noisy simulations meaned over the realizations, 1 - 
%       forward, 2 - backward) 
%   -'Q_save(2,re,nt)' (full noisy simulations for all realizations)
%   -'Qdet(2,nt)' (deterministic simulations)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
%close all; clear all;
%% 
%Integrate the 3-box model forward and backward with freshwater hosing;

% %uncomment if running here, comment if calling from Figure script
% %NOTE: usn hosing value is only valid for FAMOUS_B 1xco2
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;
tic
%load in parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
eval(params)

%select South gyre salinity exchange strength and Hosing
KS = 5.447;
H = 0;

%define timestep and gyre ramping
dt = 1*syr; %x yr in seconds
rate = .005; %Gyre ramping rate in Sv/year
Kmax = 30; % initial gyre value in Sv
Kmin = 0; % ending gyre values

tend = (Kmax-Kmin)/rate; %runtime in years
runtime = tend*syr; %seconds
nt = round(runtime/dt);
t = linspace(0,tend,nt); %years
time = linspace(0,tend/3,nt/3); %for plot (only plotting 1/3rd of total integration period)

kn1 = linspace(Kmin,Kmax,nt); %forward ramping deterministic
kn2 = linspace(Kmax,Kmin,nt); %backward ramping deterministic

%add noise
sig1 = 5e4; %set noise level
sig = sig1/sqrt(dt);  %normalize by timestep
re = 1000; %set # of realizations

N1 = sig*randn(re,nt); %Gaussian white noise with var sig^2
N2 = sig*randn(re,nt);
KN1 = repmat(kn1,re,1)+N1; %forward ramping noisy
KN2 = repmat(kn2,re,1)+N2;  %backward ramping noisy

Qusn = 0; % define shut down

%get initial conditions
solve_initial_salinities_KN

%define/save variables
KN_save(1,:) = kn1;
KN_save(2,:) = kn2;
Q_save = nan(2,re,nt);
Qvec = nan(2,nt);
Qdet = nan(2,nt);


for k = 1:2 %forward and backward
    Qvec_re = nan(re,nt);
    %set gyre strength
    if k == 1
        kn = kn1;
        KKN = KN1;
    else
        kn = kn2;
        KKN = KN2;
    end
    for p =1:re+1 %step through realizations and last run is deterministic
        %load initial salinities, resets every realization
        SN = SNsol0(k);
        ST = STsol0(k);
        
        for i = 1:nt %timestepping
            if p==re+1 %last run uses deterministic ramping KNbar
                KN = kn(i);  %deterministic ramping
            else
                KN = KKN(p,i); %stochastic ramping
            end
            
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
    
            if p == re+1
                Qdet(k,i) = Q; %save det run
            else
                Qvec_re(p,i) = Q; %save stochastic runs
            end
        end 
        
    end
    
    Qvec(k,:) = mean(Qvec_re); %avg stochastic run
    Q_save(k,:,:) = Qvec_re; %save all runs
    
end
toc
%%
%save output
save(stoch_name,'kn1','kn2','KN_save','Qvec','Q_save','Qdet','time')