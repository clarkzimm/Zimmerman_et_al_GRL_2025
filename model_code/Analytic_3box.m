% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the analytic equilibrium solutions of AMOC strength Q as a 
% function of hosing H for weak or strong gyres in the 3-box model with
% choice of model calibration. Called by Fig2_stoch_3box_stats.m,
% SI_FigS2_ratedep.m and SI_FigS3_noisedep.m
%
% Dependencies:
%   - parameters_3box_YY_ZCO2.m (parameter values for the 3-box model as 
%       calibrated to FAMOUS_B OR HadGEM2-AO (YY = [FMSB HGEM]) with Z = 
%       [1 2]xCO2) default is parameters_3box_FMSB_1CO2.m)
% 
% Output:
%   -'h1' (hosing range)
%   -'SNsol1','SNsol2' (Solution for salinity in North Atlantic box (SN) 1 -
%       positive stable, 2 - negative stable
%   -'usn','hopf','lsn' (hosing values for upper saddle node, lower saddle node, and hopf
%       bifurcations)
%   -'Qsol1','Qsol2','Qsol3' (solution for AMOC strength (Q) 1- positive
%       stable, 2- positive unstable, 3- negative stable)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024

%%
% Solve for Q analytically
%set symbolic variables
syms q gam KN KS S0 ST SN SS SIP FN FT L b H

% %uncomment if running here, comment if calling from other script
% %NOTE: usn, lsn, and hopf hosing values are only valid for FAMOUS_B 1xco2
% %choose model ('FMSB' --> FamousB, 'HGEM' --> HadGEM2-AO)
% model = 'FMSB';
% %choose CO2 level (1 --> 1xCO2 (PI-control), 2 --> 2xCO2(GW))
% CO2 = 1;
% 
%load in parameters for callibration defined above or in
%Fig2_stoch_3box_stats.m ,
% SI_FigS2_ratedep.m or SI_FigS3_noisedep.m
params = sprintf('parameters_3box_%s_%dCO2',model,CO2);
eval(params)

% % set gyre strength (1 = weak -> KN = 5.456 Sv, KS = 5.447 Sv; 2 = strong -> KN,KS = 27 Sv)
% gyre_strength = 2;

if gyre_strength == 1
    KN = 5.456; %Sv
    KS = 5.447;
    SNusn = 0.03446; %SN at the usn
    usn = 0.2138; % Hosing level at the upper saddle node bifurcation
    lsn = -0.05445; %Hosing level at the lower saddle node bifurcation
    hopf = 0.2133; %hosing level at the hopf bifurcation
    gyre = 'weak'; %for file name below
else
    if gyre_strength == 2
        KN = 27;%Sv
        KS = 27;
        SNusn = 0.03446; %SN when Q --> 0
        usn = 0.2282;  % Hosing level when Q --> 0
        gyre = 'strong'; %for file name below
    end
end

%hosing range over which solution will be calculated
h1 = linspace(-0.2,0.3,1000); 
nh = length(h1);

SNsol1 = nan(1,nh);SNsol2b=SNsol1;
SNsol2 = SNsol1;SNsol1b =SNsol1;

tic
for i = 1:length(h1) %step through each hosing value
    
    H = h1(i);

    %Q>0 equations
    eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
    eqn2 = q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt equation
    eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
    eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS)); %AMOC strength equation

    NY = vpasolve([eqn1,eqn2,eqn3,eqn4],[ST,SN,SIP,q]); %Q>0 solutions
    SNsol1(i) = max(NY.SN(1:3)); %select stable Q>0 SN solution
    SNsol2(i) = min(NY.SN(find(NY.SN>.034))); %select unstable Q>0 SN solution
    
    %Q<0 equations
    eqn1b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
    eqn2b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SN-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0;  %dST/dt equation
    eqn3b = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation

    NYb = vpasolve([eqn1b,eqn2b,eqn3b],[ST,SIP,SN]); %Q<0 solutions
    SNsol1b(i) = max(NYb.SN); %select stable Q<0 SN solution

end
toc
Qsol1 = ls.*(a*(TS-T0)+b.*(SNsol1-SS)); %Q>0 stable solution
Qsol2 = ls.*(a*(TS-T0)+b.*(SNsol2-SS)); %Q>0 unstable solution
Qsol3 = ls.*(a*(TS-T0)+b.*(SNsol1b-SS)); %Q<0 stable solution

%save the analytic solutions
fname = sprintf('analytic_sol_3box_%s_%dCO2_%s',model,CO2,gyre);
if gyre_strength == 1
    save(fname,'h1','SNsol1','SNsol2','usn','hopf','lsn','Qsol1','Qsol2','Qsol3')
else
    if gyre_strength == 2
        save(fname,'h1','SNsol1','usn','Qsol1','Qsol3') %no unstable solution for strong gyre scenario
    end
end