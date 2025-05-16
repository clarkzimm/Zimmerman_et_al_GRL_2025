% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the steady state salinities (ST,SN) for initial hosing value 
% (hmin) and final hosing value (hmax) given set parameterization/gyres. 
% Does not run independently. Called by Stochastic_sims_2box.m to provide
% initial conditions.
%
% Output:
%   -SNsol0(2,1) (initial salinity for the North Atlantic box (SN). first  
%       value for hmin, second value for hmax) 
%   -STsol0(2,1) (initial salinity for the Atlantic Thermocline box (ST). 
%       first value for hmin, second value for hmax) 
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------

%Analytic solution for set gyre strength and hosing range

SNsol0 = nan(2,1); 
for j = 1:2%find solution for hmin and hmax

    syms q ST SN  %set symbolic variables
    
    if j == 1 %solving for Q>0 (before upper saddle node)
        H = hmin;
        %Q>0 equations
        eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0;%dSN/dt equation
        eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
        eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS));%AMOC strength equation
        
        NY = vpasolve([eqn1,eqn3,eqn4],[ST,SN,q]);%Q>0 solutions
        
        SNsol0(1) = double(max(real(NY.SN(1:2))));%select stable Q>0 SN solution
    else %solving for Q<0 (after upper saddle node)
        H = hmax;
        %Q<0 equations
        eqn1b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0;%dSN/dt equation
        eqn3b = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
    
        NYb = vpasolve([eqn1b,eqn3b],[ST,SN]); %Q<0 solutions

        SNsol0(2) = double(NYb.SN(1)); %select stable Q<0 SN solution
    end 
end
