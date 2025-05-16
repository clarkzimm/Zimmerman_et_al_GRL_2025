% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the analytic equilibrium solutions of AMOC strength Q in the 
% three box model for varying North Atlantic gyre transport strength (KN) 
% with choice of model calibration. Called by SI_FigS5_rampKN.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
% 
% Output:
%   -'kn1' (range of North Atlantic gyre transport strengths)
%   -'SNsol1' (positive stable solution for salinity in North Atlantic box 
%       (SN))
%   -'Qsol1' (positive stable solution for AMOC strength (Q))
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024

% Solve for Q analytically
syms q gam KN KS S0 ST SN SS SIP FN FT L b H

% %uncomment if running here, comment if calling from other script
% %NOTE: usn, lsn, and hopf hosing values are only valid for FAMOUS_B 1xco2
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM2-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;
% 
%load in parameters for callibration defined above or in SI_FigS5_rampKN.m
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
eval(params)

%select hosing strength
H=0;
%set South Atlantic gyre strength
KS = 5.447;

%range of North ATlantic gyre strengths
kn1 = linspace(0,30,1000);
nk = length(kn1);

SNsol1 = nan(1,nk);SNsol2b=SNsol1;
SNsol2 = SNsol1;SNsol1b =SNsol1;

tic
for i = 1:length(kn1)%step through each KN
    
    KN = kn1(i);
    %Q>0 equations
    eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt eqn
    eqn2 = q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt eqn
    eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salt conservation
    eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS)); %AMOC strength

    NY = vpasolve([eqn1,eqn2,eqn3,eqn4],[ST,SN,SIP,q]); %postive Q equilibrium solutions
    SNsol1(i) = max(NY.SN); %select stable SN solution

end
toc
Qsol1 = ls.*(a*(TS-T0)+b.*(SNsol1-SS)); %solve for stable Q solution

%%
%save output
save(analytic_name,'kn1','SNsol1','Qsol1')