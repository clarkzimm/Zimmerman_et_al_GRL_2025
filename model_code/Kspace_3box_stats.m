% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the slopes of variance and autocorrelation with increasing 
% hosing for varying gyre strengths (KN,KS) in the 3-box model, calculates 
% normalized slopes as well. Can define model calibration of choice. Called
% by Fig3_Kspace_DH_stats.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - IC_Kdep_3box.m (computes the steady state salinities (SN,ST) for 
%       initial hosing value (hmin) and varying gyre strengths (KN,KS))
%   - Qdet_Kdep_3box.m (computes the deterministic trajectory of AMOC
%       strength in the 3-box model for varying gyre strengths (KN,KS).
%       Vector is cut at an interval (int) set below for lagged 
%       autocorrelation)
%          
% Output:
%   - var_slope(KN,KS) (rise in variance field as a function of KN, KS 
%       vectors)
%   - AC_slope(KN,KS) (rise in autocorrelation field as a 
%       function of KN, KS vectors)
%   - var_norm(KN,KS) (rise in variance normalized by initial variance
%       field as a function of KN, KS vectors)
%   - AC_norm(KN,KS) (rise in autocorrelation normalized by initial 
%       autocorrelation field as a function of KN, KS vectors)
% 
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------
% find the rise in statistical indicators for stochastic runs over a range
% of gyre values
%% set parameterization, timestepping, and noise level

% % uncomment if running here, comment if calling from other script
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;
% %choose lag
% int = 50; 

%load in parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
run(params)

%set timestepping
dt = 1*syr; %x yr in seconds
rate = .00005; %hosing rate in Sv/year

%set hosing range
hmax = 1; % if changed from 1 remake Qdet file (run Qdet_Kdep_3box)
hmin = -2.5; %if changed from -2.5 remake IC and Qdet file (run IC_Kdep_3box and Qdet_Kdep_3box)

%calculate timespan
tend = (hmax-hmin)/rate; %years
runtime = tend*syr; %seconds
nt = round(runtime/dt);

%set noise level
sig1 = 5e3;
sig = sig1/sqrt(dt); %normalize by timestep

h1 = linspace(hmin,hmax,nt);%deterministic hosing vector (Hbar)

re = 1000; %choose # of realizations

N1 = sig*randn(re,nt);%noise vector
H1 = repmat(h1,re,1)+N1; %stochastic hosing vector (H)

density = 100; %set # of Gyre strengths tested

%create gyre strength vectors
KNramp = linspace(2,50,density);
KSramp = linspace(2,50,density);
%create empty vectors
var_slope = nan(density,density); var_norm = var_slope;
AC_slope = nan(density,density); AC_norm = AC_slope;

% load in or solve for initial conditions for the 3-box model with model 
% calibration selected above (or in Fig3_Kspace_DH_stats.m) at hosing 
% H=hmin defined above 'IC_Kdep_3box_FMSB_1xCO2_dens100.mat' for hmin= -2.5 
% available in Zenodo in 'Output' folder
init_cond_name = sprintf('IC_Kdep_3box_%s_%dxCO2_dens%d.mat',model,CO2,density); 
if exist(init_cond_name,'file')==0
    IC_Kdep_3box %If output does not exist, takes ~25 min to run on a standard laptop Nov'24
end
load(init_cond_name)

% load in or solve for deterministic model simulations for the 3-box model 
% with model calibration selected above (or in Fig3_Kspace_DH_stats.m) for
% hosing range 'h1' and lag 'int' defined above 
% 'Qdet_Kdep_3box_FMSB_1xCO2_dens100_lag50.mat' for h1=(-2.5-->1) available 
% in Zenodo in 'Output' folder
Det_name = sprintf('Qdet_Kdep_3box_%s_%dxCO2_dens%d_lag%d.mat',model,CO2,density,int);
if exist(Det_name,'file')==0
    Qdet_Kdep_3box %If output does not exist, takes ~1 min to run on a standard laptop Nov'24
end
load(Det_name)

Qusn = 0; %define point at which near equilibrium condition no longer holds. Q=0 is when model ODEs switch

tic
for j = 1:density %step through KN values
    KN = KNramp(j);
    for l = 1:density %step through KS values
        KS = KSramp(l);
        %set initial steady-state salinities for current (KN,KS)
        SN0 = SNsol0(j,l);
        ST0 = STsol0(j,l);
        Qdet = squeeze(Qdet0(j,l,:)); %select Qdet for current (KN,KS)
        Qvec_re = nan(re,nt); %create empty stochastic Q matrix 
        for p =1:re %step through realizations
            SN = SN0;
            ST = ST0;
            for i = 1:nt %time stepping
                H = H1(p,i);
                SIP = 1/VIP*(C-(VN*SN+VT*ST+VS*SS+VD*SD)); %salinity conservation
                Q = ls*(a*(TS-T0)+b*(SN-SS)); 
                Qvec_re(p,i) = Q; 
                %update salinities
                if Q>0
                    SN = SN + (1/VN*(Q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0))*dt; 
                    ST = ST + (1/VT*(Q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0))*dt;
                elseif Q<0
                    SN = SN + (1/VN*(abs(Q)*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0))*dt;
                    ST = ST + (1/VT*(abs(Q)*(SN-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0))*dt;
                end
            end 
        end
        %compute statistics
        Qvec = mean(Qvec_re); %take average over realizations
        if mean(Qvec(:)<Qusn)>0
            L = find(Qvec(:)<Qusn,1,"first"); %find index of first shut down
        else
            L = nt; %if Q never crosses 0 define as end of vector
        end
        w = round(.05*length(1:int:L+int)); %half the window
        q_anom = Qvec_re(:,1:int:end)- Qdet'; %get anomaly
        h3 = h1(1:int:L+int); %hosing with lag accounted for
        H_ind = h3(w:end-w); %Hbar for "run"
        nH = length(H_ind);
        
        varian = nan(1,nH); rho = varian;
        for k = 1+w:nH+w-2
            q1 = q_anom(:,k+(-w:w)); %selects each window 
            q2 = q_anom(:,k+(-w:w)+1); %window one lag behind q1
            varian(k-w) = var(q2(:)); %compute variance for each window
            p = corrcoef(q2(:),q1(:)); %compute correlation between each lagged window
            rho(k-w) = p(1,2); %select autocorrelation
        end
        %compute slope
        run_var = H_ind(nH - 3*w)-H_ind(round(nH/5)); %hosing range over which variance slope is calculated
        run_AC = H_ind(end - 2*w)-H_ind(1); %hosing range over which AC slope is calculated
        rise_var = varian(nH - 3*w)-varian(round(nH/5)); %rise in variance
        rise_AC = rho(end - 2*w)-rho(1); %rise in AC
        %slopes
        var_slope(j,l) = rise_var/run_var; 
        AC_slope(j,l) = rise_AC/run_AC;
        %normalized slopes
        var_norm(j,l) = var_slope(j,l)/varian(round(nH/5));
        AC_norm(j,l) = AC_slope(j,l)/rho(1);  
    end
end
toc
%%
fname = sprintf('stats_Kdep_3box_%s_%dCO2_lag%d.mat',model,CO2,int); 
save(fname,'AC_slope','var_slope','var_norm','AC_norm','KNramp','KSramp')
