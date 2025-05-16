% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Computes hysteresis width (Delta H*) in hosing (H) for deterministic runs 
% of the 3-box model for varying gyre strengths (KN,KS). Called by
% SI_FigS6_dethyst.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
%   - IC_Kdep_3box.m (computes the steady state salinities (SN,ST) for 
%       initial hosing value (hmin) and varying gyre strengths (KN,KS))
%          
% Output:
%   - detDH(KN,KS) (hysteresis width field as functions of KN, KS vectors)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------

%% 
% % uncomment if running here, comment if calling from other script
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;

%load parameter values
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
run(params)

%define timestepping and hosing rate
dt = 1*syr; %x yr in seconds
rate = .00005; %Sv/year

%set hosing range
hmax = 1;
hmin = -2.5;

tend = (hmax-hmin)/rate; %years
runtime = tend*syr; %seconds
nt = round(runtime/dt);

%deterministic hosing vectors 
h1 = linspace(hmin,hmax,nt); %forward
h2 = linspace(hmax,hmin,nt); %backward

%set empty arrays with data density of choice
density = 100;
KNramp = linspace(2,50,density);
KSramp = linspace(2,50,density);
detDH = nan(density,density);
SNsol0 = nan(density,density);
STsol0 = nan(density,density);

% load in or solve for initial conditions for the 3-box model with model 
% calibration selected above (or in SI_FigS6_dethyst.m) at hosing 
% H=hmin defined above 'IC_Kdep_3box_FMSB_1xCO2_dens100.mat' for hmin= -2.5 
% available in Zenodo in 'Output' folder
init_cond_name = sprintf('IC_Kdep_3box_%s_%dxCO2_dens%d.mat',model,CO2,density); 
if exist(init_cond_name,'file')==0
    IC_Kdep_3box %If output does not exist, takes ~25 min to run on a standard laptop Nov'24
end
load(init_cond_name)

Qusn = 3; %set where to measure hysteresis
tic
for j = 1:density %step through KN values
for l = 1:density %step through KS values
    KN =KNramp(j);
    KS =KSramp(l);
    %load initial conditions
    SN0 = SNsol0(j,l);
    ST0 = STsol0(j,l);
    H_bif = nan(2,1);

    for k = 1:2 %run forward and backward
        Qvec = nan(nt,1);
        if k==1 %forward
            SN = SN0; %equilibrium salinities
            ST = ST0;
            h = h1;
        else %backward
            SN = SN_end; %final salinities from forward run
            ST = ST_end;
            h = h2;
        end
        for i = 1:nt %integrate
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
            Qvec(i) = Q;
        end
        SN_end = SN; %final salinities approximately steady-state IC for reverse hosing
        ST_end = ST;
      
        if k ==1 %forward
            if mean(Qvec(:)<Qusn)>0 %ensure system does shut down
                H_bif(k) = h(find(Qvec(:)<Qusn,1,"first")); %find hosing when simulation crosses Q=3
            end
        else
            if mean(Qvec(:)>Qusn)>0 %ensure system does turn on
                H_bif(k) = h(find(Qvec(:)>Qusn,1,"first")); %find hosing when simulation crosses Q=3
            end
        end
    end    
    if anynan(H_bif)==0 %check for any which did not complete a hysteresis loop
        detDH(j,l) = H_bif(1)-H_bif(2);%calculate hysteresis width
    else
        detDH(j,l) = 10; %in case the system does not jump to lower (or upper) branch set hysteresis value to 10 (much larger than actual hysteresis width)
        fprintf('(%d,%d);',j,l)
    end
end
end
toc
%%
%file name defined in SI_FigS6_dethyst.m
save(det_name,'detDH','KNramp','KSramp')