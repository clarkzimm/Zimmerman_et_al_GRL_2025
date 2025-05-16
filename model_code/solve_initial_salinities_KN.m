% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the steady state salinities (ST,SN) for initial gyre value 
% (knmax) and final gyre value (knmin) given set parameterization/hosing. 
% Does not run independently. Called by Stochastic_sims_3box_rampKN.m to 
% provide initial conditions.
%
% Output:
%   -SNsol0(2,1) (initial salinity for the North Atlantic box (SN). first  
%       value for knmin, second value for knmax) 
%   -STsol0(2,1) (initial salinity for the Atlantic Thermocline box (ST). 
%       first value for knmin, second value for knmax) 
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------

%steady state analytic solution for set hosing strength and gyre range

STsol0 = nan(2,1); SNsol0 = STsol0;
for j = 1:2
    syms q ST SN SIP
    if j == 1
        KN = Kmin;
    else
        KN = Kmax;
    end 
    %Q>0 equations: AMOC is always "on"
    eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
    eqn2 = q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt equation
    eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
    eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS));%AMOC strength equation
    
    NY = vpasolve([eqn1,eqn2,eqn3,eqn4],[ST,SN,SIP,q]);%Q>0 solutions
    
    STsol0(j) = min(real(NY.ST(find(NY.ST>.0354)))); %select stable Q>0 ST solution
    SNsol0(j) = max(real(NY.SN(1:3))); %select stable Q>0 SN solution
end

