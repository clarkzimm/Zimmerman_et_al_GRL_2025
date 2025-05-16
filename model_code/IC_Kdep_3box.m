% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Computes the steady state salinities (SN,ST) for initial hosing value 
% (hmin) and varying gyre strengths (KN,KS). Does not run independently. 
% Called by Stochastic_3box_DHstar.m, Kspace_3box_stats.m, and 
% Det_DH_3box.m to provide initial conditions.
%          
% Output:
%   - SNsol0(KN,KS) (initial North Atlantic salinity (SN) field as a 
%       function of KN, KS vectors)
%   - STsol0(KN,KS) (initial Atlantic Thermocline salinity (ST) field as a 
%       function of KN, KS vectors)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------

%Analytic solution for initial salinities across a range of Kn and Ks values

SNsol0 = nan(density,density);
STsol0 = nan(density,density);
tic
for j = 1:density
for l = 1:density

    KN =KNramp(j); %set gyre strengths; KNramp and KSramp vectors defined in Stochastic_3box_DHstar.m and Kspace_3box_stats.m
    KS =KSramp(l);

    syms q ST SN SIP %set symbolic variables
    
    H = hmin; %initial hosing value
    
    %Q>0 equations: always starting from "on" AMOC
    eqn1 = q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*H)*S0 == 0; %dSN/dt equation
    eqn2 = q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*H)*S0 == 0; %dST/dt equation
    eqn3 = C == VN*SN+VT*ST+VS*SS+VIP*SIP+VD*SD; %salinity conservation equation
    eqn4 = q == ls*(a*(TS-T0)+b*(SN-SS)); %AMOC strength equation
    
    NY = vpasolve([eqn1,eqn2,eqn3,eqn4],[ST,SN,SIP,q]); %Q>0 solutions
    
    STsol0(j,l) = double(max(real(NY.ST(find(NY.ST<.039))))); %select stable Q>0 ST solution; convert from symbolic variable to standard variable
    SNsol0(j,l) = double(max(NY.SN(1:3))); %select stable Q>0 SN solution; convert from symbolic variable to standard variable

end
end
toc
%init_cond_name defined in Stochastic_3box_DHstar.m and Kspace_3box_stats.m
save(init_cond_name,'STsol0','SNsol0') %save matrix of initial conditions
