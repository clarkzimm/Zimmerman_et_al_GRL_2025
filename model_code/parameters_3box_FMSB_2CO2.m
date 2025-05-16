% Code accompanying Zimmerman et al (GRL, 2024)
%
% Summary: 
% Parameter values for the 3-box model as calibrated to FAMOUS_B 2xCO2
% -------------------------------------------------------------------------
% Clark Zimmerman (cczimmerman3@wisc.edu) 
% Till Wagner (till.wagner@wisc.edu)
% November 2024
% -------------------------------------------------------------------------
syr = 31536926; %s/yr convert between seconds and years

%Volumes of the boxes
VN = 3.683*1e10; %Sv/s %note this is FAMOUS_B
VT = 5.418*1e10;
VS = 6.097*1e10;
VIP = 14.860*1e10;
VD = 99.25*1e10;

%intitial salinities
S0 = .035; %psu  reference salinity
SS = 0.034427; % fixed in the 3 box model
SD = 0.034538; % fixed in the 3 box model
SN0 = 0.034912;
ST0 = 0.035435;
SIP0 = 0.034668;
C = VN*SN0+VT*ST0+VS*SS+VIP*SIP0+VD*SD;

%air-sea freshwater flux
FN0 = 0.486; %Sv
FT0 = -0.997;

%Temperatures
TS = 7.919; %C
T0 = 3.87;

%constants
mu = 22 * 1e-2;% Cm^-3s
l = 1.62 * 1e1;% Sv m^3kg^-1
gam = 0.36;%unitless
a = 0.12; %kgm^-3K^-1
b = 790; %kg m^-3 psu^-1
ls = l/(1+l*a*mu);

%hosing fraction
FNhos = .1311;
FThos = .6961;

