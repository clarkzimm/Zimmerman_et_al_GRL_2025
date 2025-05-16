% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Parameter values for the 5-box model as calibrated to AOGCM HadGEM2-AO 2xCO2
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
syr = 31536926; %s/yr convert between seconds and years

%Volumes of the boxes
VN = 5.259*1e10; %Sv/s
VT = 7.4*1e10;
VS = 9.336*1e10;
VIP = 19.220*1e10;
VD = 89.90*1e10;

%intitial salinities
S0 = .035; %psu  reference salinity
SS = 0.034427; % fixed in the 3 box model
SD = 0.034538; % fixed in the 3 box model
SN0 = 0.0351036;
ST0 = 0.035642;
SIP0 = 0.0345665;
C = VN*SN0+VT*ST0+VS*SS+VIP*SIP0+VD*SD;

%air-sea freshwater flux
FN0 = 0.496; %Sv
FT0 = -0.921;
FS = 1.021;
FIP = -0.596;

%Temperatures
TS = 7.424; %C
T0 = 3.29;

%constants
mu = 16 * 1e-2;% Cm^-3s
l = 1.66 * 1e1;%m^6kg^-1s^-1
gam = 0.73;%unitless
eta = 9.871; %Sv
a = 0.12; %kgm^-3K^-1
b = 790; %kg m^-3 psu^-1
ls = l/(1+l*a*mu);

%hosing fraction
hosF_N = .285;
hosF_T = .522;
hosF_S = -0.299;
hosF_IP = -0.508;


KIP = 1029.641; %Sv


