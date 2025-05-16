% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Solves for the steady state salinities (ST,SN) for initial hosing value 
% (hmin) and final hosing value (hmax) given set parameterization/gyres. 
% Does not run independently. Called by Stochastic_sims_3box.m  and 
% SI_FigS2_ratedep.m to provide initial conditions.
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

%steady state analytic solution for set gyre strength and hosing values

STsol0 = nan(2,1); SNsol0 = STsol0;
H0 = [hmin hmax];

for j = 1:2 %find solution for hmin and hmax

    syms q ST SN SIP %set symbolic variables
    H = H0(j);
    if H < usn %determine if solving for Q>0 (before upper saddle node)
        %Q>0 equations
        eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
        eqn2 = q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt equation
        eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
        eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS));%AMOC strength equation
        
        NY = vpasolve([eqn1,eqn2,eqn3,eqn4],[ST,SN,SIP,q]);%Q>0 solutions
        
        STsol0(j) = min(real(NY.ST(find(NY.ST>.0354)))); %select stable Q>0 ST solution
        SNsol0(j) = max(real(NY.SN(1:3))); %select stable Q>0 SN solution
    else %solving for Q<0 (after upper saddle node)
        %Q<0 equations
        eqn1b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
        eqn2b = abs(ls*(a*(TS-T0)+b*(SN-SS)))*(SN-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt equation
        eqn3b = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
    
        NYb = vpasolve([eqn1b,eqn2b,eqn3b],[ST,SIP,SN]); %Q<0 solutions
        SNsol0(j) = NYb.SN(1); %select stable Q<0 SN solution
        STsol0(j) = NYb.ST(1); %select stable Q<0 ST solution
    end 
end

