% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Parameter values for the 5-box model as calibrated to FAMOUS_B 1xCO2
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
syr = 31536926; %s/yr convert between seconds and years

%Volumes of the boxes
VN = 3.261*1e10; %Sv*s %note this is FAMOUS_B
VT = 7.777*1e10;
VS = 8.897*1e10;
VIP = 22.020*1e10;
VD = 86.490*1e10;

%intitial salinities
S0 = .035; %psu  reference salinity
SS = 0.034427;
SD = 0.034538; 
SN0 = 0.034912;
ST0 = 0.035435;
SIP0 = 0.034668;
C = VN*SN0+VT*ST0+VS*SS+VIP*SIP0+VD*SD;

%air-sea freshwater flux
FN0 = 0.384; %Sv
FT0 = -0.723;
FS = 1.078;
FIP = -0.739;

%Temperatures
TS = 4.773; %C
T0 = 2.650;

%constants
mu = 5.5 * 1e-2;% Cm^-3s
l = 2.79 * 1e1;% Sv m^3kg^-1
gam = 0.39;%unitless
eta = 74.492; %Sv
a = 0.12; %kgm^-3K^-1
b = 790; %kg m^-3 psu^-1
ls = l/(1+l*a*mu);

%hosing fraction
hosF_N = .070;
hosF_T = .752;
hosF_S = -0.257;
hosF_IP = -0.565;

KIP = 96.817;



