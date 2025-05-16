% Code accompanying Zimmerman et al (GRL, 2025)
%
% Summary: 
% Solves for the analytic equilibrium solutions of AMOC strength Q as a 
% function of hosing H for weak or strong gyres in the 2-box model with 
% choice of model calibration. Called by SI_FigS4_2box.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
% 
% Output:
%   -'han' (hosing range)
%   -'usn','lsn' (hosing values for upper saddle node and lower saddle node
%       bifurcations)
%   -'Qsol1','Qsol2','Qsol3' (solution for AMOC strength (Q) 1- positive
%       stable, 2- positive unstable, 3- negative stable)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% January 2025

%%
tic
% Solve for Q analytically
syms q ST SN 

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
    usn = 0.3323; % Hosing level at the upper saddle node bifurcation
    lsn = -1.7688; %Hosing level at the lower saddle node bifurcation
    han = linspace(-2,1,1000); %set hosing range
else
    if gyre_strength == 2
        KN = 22;%Sv
        usn = 9.5018;  % Hosing level when Q --> 0
        han = linspace(7,10,1000); %set hosing range
    end
end

nh = length(han);

SNsol1 = nan(1,nh);
SNsol2 = SNsol1;SNsol1b =SNsol1;

for i = 1:length(han)  %step through each hosing value
    
    H = han(i);
    %Q>0 equations
    eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0;
    eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD;%salinity conservation equation
    eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS));%AMOC strength equation

    NY = vpasolve([eqn1,eqn3,eqn4],[ST,SN,q]);%Q>0 solutions
    SNsol1(i) = max(NY.SN(:)); %select stable Q>0 SN solution
    SNsol2(i) = min(NY.SN(:)); %select unstable Q>0 SN solution
    
    %Q<0 equations
    eqn1b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0;%dSN/dt equation
    eqn3b = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD;%salinity conservation equation

    NYb = vpasolve([eqn1b,eqn3b],[ST,SN]);%Q<0 solutions
    SNsol1b(i) = max(NYb.SN);%select stable Q<0 SN solution

end

Qsol1 = ls.*(a*(TS-T0)+b.*(SNsol1-SS)); %Q>0 stable solution
Qsol2 = ls.*(a*(TS-T0)+b.*(SNsol2-SS)); %Q>0 unstable solution
Qsol3 = ls.*(a*(TS-T0)+b.*(SNsol1b-SS)); %Q<0 stable solution

%save the analytic solutions
if gyre_strength == 1
    save(analytic_name,'han','usn','lsn','Qsol1','Qsol2','Qsol3')
else
    if gyre_strength == 2
        save(analytic_name,'han','usn','Qsol1','Qsol3') %no unstable solution for strong gyre scenario
    end
end
toc