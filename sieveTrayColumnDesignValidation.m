%Sieve Tray Column Script for Acetone n-butanol ethanol system
Spacing = [9 12 18 24 36];

%calculated with an overall mass balance mass balance previously. Densities
%and average MW through gas laws or literature values for the system at set
%mole fractions

V=[674 162.96] %top, bottom
L =[285 774]
rhoL = [736.3940979 758.3217177]
rhoV = [2.091852323 2.193733207]
MWL = [63.05271 72.67731]
MWV = [60.64656 63.60027]
 
F = ((L.*MWL)/(V.*MWV))*sqrt(rhoV./rhoL);
Csbf = [(10.^(-1.0674-0.55780.*log10(F)-0.17919.*(log10(F)).^2)) (10.^(-1.1622-0.56014.*log10(F)-0.18168.*(log10(F)).^2)) (10.^(-1.0262-0.63513.*log10(F)-0.20097.*(log10(F)).^2)) (10.^(-0.94506-0.70234.*log10(F)-0.22618.*(log10(F)).^2))*(10.^(-0.85984-0.73980.*log10(F)-0.23735)).*(log10(F)).^2]
surfaceTensionMix =17.3 %literature value for the system
 
uFloodTOP = (Csbf.*(surfTenMix/20).^0.2).
*sqrt((rhoL(1)-rhoV(1))./(rhoV(1)))
DTOP = sqrt((4*388.889*8.314*(80.15+273.15))./
(3600.*3.1415.*0.9.*0.75.*uFloodTOP));
uFloodBOT = (Csbf.*(surfTenMix/20).^0.2).
*sqrt((rhoL(2)-rhoV(2))./(rhoV(2)))
DBOT = sqrt((4*388.889*8.314*(80.15+273.15))./
(3600.*3.1415.*0.9.*0.75.*uFloodBOT));
Dratio = DTOP./DBOT;
 
HeightCol = (6.*Spacing)./12
%Hydraulics. assume 75% flooding. Plate top and bottom within 1% same, calculate with largest F
Entrainment = 0.09 %from plot 
absoluteEntrainment = (Entrainment.*L(1))./(1-Entrainment)
TotalFlow = (absoluteEntrainment+L(1)) %These values of entrainment are reasonable.
the typical Iweir/Dia falls in range 0.6-0.75,
 take it to be 0.726 giving n to be 0.9
A_active =(pi().*DTOP.^2./4)*(2*0.9-1); 
%beta to be 0.1
A_hole = 0.1*A_active
nHoles = ceil(A_hole./(pi().*((0.5/12)^2/4)))
VapourVelHole=(V(1).*MWV(1))./(3600.*(rhoV(1)./16.01).*(A_hole))
Lg=TotalFlow*MWL(1)*(1/41.12)*7.48*(1/60);

 %Orifice coeff. Take a standard 14 gauge tray, 
this do/t_tray = 2.4.
Co=0.85032-0.04231*2.4+0.0017954*2.4^2;
hdp_dry = 0.003.*(VapourVelHole.^2).*(rhoV(1)./16.01).*((61.03./rhoV(1)).*(1-0.1^2)/Co.^2)

%assume 1 in gap. Adu = (1/12).*(0.726.*DTOP);
Abscissa = Lg./(0.726.*DTOP).^2.5; %weir correction is 1 by graphical. 
h_crest = 0.092.*(Lg./Adu).^(2/3);%Assume h_grad is zero, H_weir is 2 inches
hdu = 0.56.*(Lg./(449.*Adu)).^2;
hdc=hdu+hdp_dry+h_crest+2;
hdcAerated = hdc./0.5; %hdc Aerated is half the spacing for 12 in. 
Ad = (pi().*DTOP.^2./4).*0.1;
t_res= (Ad.*hdc.*3600.*(rhoL(1)./16.0185))./(TotalFlow.*MWL(1)*12); 
%Greater than 3s, residence time is good.
x = 2+h_crest;
INR = 0.10392+0.25119.*x-0.021675.*x.^2;
ho = (0.040*surfTenMix)/((rhoL(1)./16.0185)*0.5);
INL = hdp_dry+ho;
if INL>INR
    disp("No Weeping")
end

