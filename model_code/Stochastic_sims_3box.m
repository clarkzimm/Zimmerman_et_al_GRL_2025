% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Integrates the 3-box model for both forward and backward hosing (H),
% solving for AMOC strength (Q), with Guassian white noise added to the
% hosing. Called by Fig2_stoch_3box_stats.m and SI_FigS3_noisedep.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - solve_initial_salinities.m (solves for the steady state salinities 
%       (SN,ST) for initial hosing value (hmin) and final hosing value (hmax)
%          
% Output:
%   -'h_save','h1','h2' (deterministic hosing (bar(H)), h1- forward, h2-
%       backward, h_save = [h1 h2]
%   -'Qvec(2,nt)' (noisy simulations meaned over the realizations, 1 - forward, 2 - backward) 
%   -'Q_save(2,re,nt)' (full noisy simulations for all realizations)
%   -'Qdet(2,nt)' (deterministic simulations)
%   -'H_bif_min' (index of first crossing of Qusn, defined below (default =
%       0)
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
% 
% % set gyre strength (1 = weak -> KN = 5.456 Sv, KS = 5.447 Sv; 2 = strong -> KN,KS = 27 Sv)
% gyre_strength = 1;

%load in parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
eval(params)

if gyre_strength == 1
    KN = 5.456; %Sv
    KS = 5.447;
    SNusn = 0.03446; %SN at the usn
    usn = 0.2138; % Hosing level when Q --> 0
    gyre = 'weak';
else
    if gyre_strength == 2
        KN = 27;
        KS = 27;
        SNusn = 0.03446;
        usn = 0.2282;  % Hosing level when Q --> 0
        gyre = 'strong';
    end
end

%define timestep and hosing
dt = 1*syr; %x yr in seconds
rate = .00005; %Hosing rate in Sv/year
hmax = 1; % ending hosing values
hmin = -1; %initial hosing value

tend = (hmax-hmin)/rate; %runtime in years
runtime = tend*syr; %seconds
nt = round(runtime/dt);

h1 = linspace(hmin,hmax,nt); %forward hosing deterministic
h2 = linspace(hmax,hmin,nt); %backward hosing deterministic
h_save = nan(2,nt);

%add noise
sig1 = 5e3; %set noise level
sig = sig1/sqrt(dt); %normalize by timestep
re = 1000; %set # of realizations

N1 = sig*randn(re,nt); %Gaussian white noise with var sig^2
N2 = sig*randn(re,nt);
H1 = repmat(h1,re,1)+N1; %forward hosing noisy
H2 = repmat(h2,re,1)+N2; %backward hosing noisy

Qusn = 0; % define shut down

%get initial conditions
solve_initial_salinities

%define/save variables
h_save(1,:) = h1;
h_save(2,:) = h2;
Q_save = nan(2,re,nt);
Qvec = nan(2,nt);
Qdet = nan(2,nt);
H_bif = nan(re,1);H_bif2 = H_bif;H_bif1=H_bif;

tic
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
        SN = SNsol0(k); %load initial salinities, resets every realization
        ST = STsol0(k);

        for i = 1:nt %timestepping
            if p==re+1 %last run uses deterministic hosing Hbar
                H = h(i); %deterministic hosing
            else
                H = HH(p,i); %noisy hosing
            end
            
            SIP = 1/VIP*(C-(VN*SN+VT*ST+VS*SS+VD*SD));%salinity conservation
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

        if p < re+1 %only for stochastic runs with increasing hosing
            if k == 1
                H_bif(p) = find(Qvec_re(p,:)<Qusn,1,"first"); %index of shut down
            end
        end
    end

    H_bif_min = min(H_bif); %earliest bif across realizations
    Qvec(k,:) = mean(Qvec_re); %avg stochastic run
    Q_save(k,:,:) = Qvec_re; %save all runs

end
toc

% save output
fname = sprintf('stochastic_3box_%s_%dCO2_%s',model,CO2,gyre);
save(fname,'h1','h2','h_save','Qvec','Q_save','Qdet','H_bif_min')