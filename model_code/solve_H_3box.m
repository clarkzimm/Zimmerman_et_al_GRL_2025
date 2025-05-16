% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% analytic equilibrium solutions of hosing (H) as a function of AMOC 
% transport (q), and gyre transport strengths (KN,KS) for 3-box model. Does 
% not run independently. Called by Analytic_DeltaH.m
%
%          
% Output:
%   - Hpos(q,KN,KS) (solution of hosing (H) for positive AMOC transport (q))
%   - Hneg(q,KN,KS) (solution of hosing (H) for negative AMOC transport (q))
%
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% October 2024
% -------------------------------------------------------------------------
% solutions of H as a function of q (from Mathematica) for 3-box model

%redefine variables for simpler equations
A = ls*a; B = ls*b; b2 = B*(1-gam);d = A*(TS-T0)-B*SS; d2 = d*(1-gam);
rn = VN/VIP;rt = VT/VIP; z = (C-SD*VD-SS*VS)/VIP;
P2 = (d-B*SD)/VN; P4 = (d*SD+FN0*S0)/VN; P7 = B/VT;P9 = -FT0*S0/VT;
C2 =B/VN; C3 = -S0*FNhos/VN;C13 = -S0*FThos/VT;C4 = S0*FN0/VN; C5 = (B*gam*SS-d2*rn+b2*z)/VT;
C6 = b2*rn/VT; C7 = (b2*rt+B)/VT; C8 = (d2*rt+d)/VT; C9 = (SS*d*gam+d2*z-FT0*S0)/VT;

%steady state solution for hosing (H) when AMOC strength (Q) is positive
Hpos=@(q,KN,KS) (((-C2)*(C6+C7).*(d-q).^3+B.*(d-q).^2.*((C6+C7).*(d+KN)/VN+C2.*(-C5+C8+KS/VT))+B^3.*(C4.*(C8+(KN+KS)./VT)-...
    (d+KN)./VN.*(C9+KS.*SS./VT))+B^2.*(d-q).*((-C4)*C7+(d+KN)./VN.*(C5-C8-KS./VT)+C2.*(C9+KS.*SS./VT)))./...
  (B^2.*(B*C13.*(d+KN)./VN-(C13*C2+C3*C7).*(d-q)+B*C3.*(C8+(KN+KS)./VT))));
%steady state solution for hosing (H) when AMOC strength (Q) is negative
Hneg=@(q,KN,KS) (((-C2).*(d-q).^2.*((-P7).*q+B.*(KN+KS)./VT+d*(P7-B/VT))+B.*(P2*P7.*q.^2+d^2*P2*(P7-B/VT)-... 
      B.*q.*(P4*P7+KS.*P2./VT+KN.*((-KS)./VN+P2)./VT)+d.*(-2*P2*P7.*q-B^2*P4/VT+... 
         B.*(P4*P7+KN.*((-KS)./VN+P2)./VT+P2.*(KS+q)./VT))+B^2.*(KN.*P4./VT+KS.*P4./VT-KN./VN.*(P9+KS.*SS./VT))))./...
   (B^2.*(B*C13.*KN./VN+C3*P7.*(d-q)+B*C3.*(-d+KN+KS)./VT)));