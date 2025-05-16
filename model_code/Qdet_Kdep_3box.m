% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% computes the deterministic trajectory of AMOC strength in the 3-box model 
% for varying gyre strengths (KN,KS). Vector is cut at an interval (int) 
% set in Kspace_3box_stats.m. Does not run independently. Called by 
% Kspace_3box_stats.m
% 
%          
% Output:
%   - Qdet0(KN,KS) (deterministic AMOC strength (Q) trajectory at lag (int)
%       field as a function of KN,KS vectors)
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------
Qdet01 = nan(density,density,nt); %full deterministic trajectory vector
Qdet0 = nan(density,density,nt/int); % deterministic trajectory vector cut at lag 'int'
tic
for j = 1:density
    KN = KNramp(j);
    for l = 1:density
        KS = KSramp(l);
        SN = SNsol0(j,l); %load in initial salinities for each (KN,KS)
        ST = STsol0(j,l);
        for i = 1:nt %time stepping
            h = h1(i);
            SIP = 1/VIP*(C-(VN*SN+VT*ST+VS*SS+VD*SD)); %salt conservation
            Q = ls*(a*(TS-T0)+b*(SN-SS));
            if Q>0
                dSN = 1/VN*(Q*(ST-SN)+KN*(ST-SN)-(FN0+FNhos*h)*S0);
                dST = 1/VT*(Q*(gam*SS+(1-gam)*SIP-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*h)*S0);
            elseif Q<0
                dSN = 1/VN*(abs(Q)*(SD-SN)+KN*(ST-SN)-(FN0+FNhos*h)*S0);
                dST = 1/VT*(abs(Q)*(SN-ST)+KS*(SS-ST)+KN*(SN-ST)-(FT0+FThos*h)*S0);
            end
            SN = SN + dSN*dt; %update salinities
            ST = ST + dST*dt; 
            Qdet01(j,l,i) = Q; %fill full AMOC trajectory
        end
        Qdet0(j,l,:) = Qdet01(j,l,1:int:nt); %cut vector at specified interval
    end
end
toc
%Det_name set in Kspace_3box_stats.m
save(Det_name,'Qdet0') %save cut deterministic trajectory