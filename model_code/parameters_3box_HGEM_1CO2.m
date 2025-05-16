% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Parameter values for the 3-box model as calibrated to AOGCM HadGEM2-AO 1xCO2
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
syr = 31536926; %s/yr convert between seconds and years

%Volumes of the boxes
VN = 3.557*1e10; %Sv/s 
VT = 8.908*1e10;
VS = 10.330*1e10;
VIP = 19.219*1e10;
VD = 90.23*1e10;

%intitial salinities
S0 = .035; %psu  reference salinity
SS = 0.034427; % fixed in the 3 box model
SD = 0.034538; % fixed in the 3 box model
SN0 = 0.0351036;
ST0 = 0.035642;
SIP0 = 0.0345665;
C = VN*SN0+VT*ST0+VS*SS+VIP*SIP0+VD*SD;

%air-sea freshwater flux
FN0 = 0.453; %Sv
FT0 = -0.798;

%Temperatures
TS = 6.456; %C
T0 = 2.71;

%constants
mu = 1.4 * 1e-2;% Cm^-3s
l = 2.17 * 1e1;%m^6kg^-1s^-1
gam = 0.85;%unitless
a = 0.12; %kgm^-3K^-1
b = 790; %kg m^-3 psu^-1
ls = l/(1+l*a*mu);

%hosing fraction
FNhos = .117;
FThos = .703;

