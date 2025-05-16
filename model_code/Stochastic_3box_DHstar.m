% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Computes hysteresis width (Delta H*) in hosing (H) for stochastic runs 
% of the 3-box model for varying gyre strengths (KN,KS). Called by
% Fig3_Kspace_DH_stats.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - IC_Kdep_3box.m (computes the steady state salinities (SN,ST) for 
%       initial hosing value (hmin) and varying gyre strengths (KN,KS))
%          
% Output:
%   - DH(KN,KS) (hysteresis width field as functions of KN, KS vectors)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------

%% set parameterization

% % uncomment if running here, comment if calling from other script
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;

%load parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
run(params)
%% set timestepping and noise level
%set timestepping
dt = 1*syr; %x yr in seconds
rate = .00005; %hosing rate in Sv/year

%set hosing range
hmax = 1; % ending hosing value
hmin = -2.5; %initial hosing value; if changed from -2.5 remake IC file (run IC_Kdep_3box)

%calculate timespan
tend = (hmax-hmin)/rate; %years
runtime = tend*syr; %seconds
nt = round(runtime/dt);

%set noise level
sig1 = 5e3;
sig = sig1/sqrt(dt); %normalize by timestep

h1 = linspace(hmin,hmax,nt);%deterministic hosing vector (Hbar) for forward ramping
h2 = linspace(hmax,hmin,nt);%deterministic hosing vector (Hbar) for backward ramping

re = 1000; %choose # of realizations

%noise vectors
N1 = sig*randn(re,nt);
N2 = sig*randn(re,nt);
%stochastic hosing vector (H)
H1 = repmat(h1,re,1)+N1;
H2 = repmat(h2,re,1)+N2;

density = 100; %set # of Gyre strengths tested
%create gyre strength vectors
KNramp = linspace(2,50,density);
KSramp = linspace(2,50,density);

%create empty matrices
DH = nan(density,density);
SNsol0 = nan(density,density);
STsol0 = nan(density,density);

% load in or solve for initial conditions for the 3-box model with model 
% calibration selected above (or in Fig3_Kspace_DH_stats.m) at hosing 
% H=hmin defined above 'IC_Kdep_3box_FMSB_1xCO2_dens100.mat' for hmin= -2.5 
% available in Zenodo in 'Output' folder
init_cond_name = sprintf('IC_Kdep_3box_%s_%dxCO2_dens%d.mat',model,CO2,density); 
if exist(init_cond_name,'file')==0
    IC_Kdep_3box %If output does not exist, takes ~25 min to run on a standard laptop Nov'24
end
load(init_cond_name)

Qusn = 3; %set where to measure hysteresis

for j = 1:density %step through KN values
for l = 1:density %step through KS values
    KN =KNramp(j);
    KS =KSramp(l);
    %set initial steady-state salinities for current (KN,KS)
    SN0 = SNsol0(j,l);
    ST0 = STsol0(j,l);
    H_bif = nan(2,1);
    %create vectors for ending salinities
    SN_end = nan(re,1); 
    ST_end = nan(re,1);
tic
    for k = 1:2 %run forward then backward
        Qvec_re = nan(re,nt);
        
        for p =1:re
            if k==1 %forward
                SN = SN0;% initialize with steady state values
                ST = ST0;
                h = h1;
            else %backward
                SN = mean(SN_end); %initialize with mean ending values
                ST = mean(ST_end);
                h = h2;
            end
            
            for i = 1:nt %timestepping
                if k==1 %forward
                    H = H1(p,i);
                else %backward
                    H = H2(p,i);
                end
                
                SIP = 1/VIP*(C-(VN*SN+VT*ST+VS*SS+VD*SD)); %salinity conservation
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
        
                Qvec_re(p,i) = Q;
            end 
            %on last forward timestep fill end vectors with final
            %salinities
            if k == 1
                SN_end(p) = SN;
                ST_end(p) = ST;
            end
        end
        Qvec = mean(Qvec_re); %take average over realizations
        if mean(Qvec(:)<Qusn)>0 %check to be sure Q-->0 when forwards
            if k == 1
                H_bif(k) = h1(find(Qvec(:)<Qusn,1,"first")); %find hosing when average realization crosses Q=3
            else 
                if mean(Qvec(:)>Qusn)>0 %check to be sure Q --> + when backwards
                    H_bif(k) = h2(find(Qvec(:)<Qusn,1,"last")); %find hosing when average realization crosses Q=3
                end
            end
        end
    end
toc    
    if anynan(H_bif)==0 %check that system completed a hysteresis loop
        DH(j,l) = H_bif(1)-H_bif(2); %calculate hysteresis width
    else
        DH(j,l) = 10; %in case the system does not jump to lower (or upper) branch set hysteresis value to 10 (much larger than actual hysteresis width)
    end
end
end
%%
fname = sprintf('stoch_DH_3box_%s_%dCO2.mat',model,CO2);
save(fname,'DH','KNramp','KSramp')

