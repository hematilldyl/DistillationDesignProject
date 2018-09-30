clear
clc

%Mass balances
F = 40;
R = 2.5;
V = 30;
D = V/1.5;
B = F-D;
xM = ((20/32.04))/((20/32.04)+(80/60.1));
yM = 0.41; %from xy
alphaMI = (yM*(1-xM))/(xM*(1-yM));
slopeROL= R/(R+1);

Mcp = @(T) 40.152*T+(3.1046E-1)*T-(1.0291E-3)*T.^2-(1.4598E-6)*T.^3;
Icp = @(T) 72.525*T+(7.9553E-1)*T-(2.6330E-3)*T.^2-(3.6498E-6)*T.^3;

data = xlsread('hxydata.xlsx');

x = data(:,1);
y = data(:,2);
T = data(:,3);
figure;
%xy
plot(x,y,x,x),xlabel("x Fraction Methanol"),...
    ylabel("y Fraction n-Methanol"),title("xy Diagram n-Propanol/Isopropanol")

figure;
%txy
plot(x,T,y,T),xlabel("x/y Fraction Methanol"),...
    ylabel("Temperature (Deg C)"),title("Txy Diagram n-Propanol/Isopropanol")


HvapM = 35140;
HvapI = 39870;
%from txy
Tbp = 75.8+273.15;
Tfeed = 273.15+20;
vapMixFeed = xM*HvapM+(1-xM)*HvapI;
q = ((xM.*integral(Mcp,Tfeed,Tbp)+((1-xM).*integral(Icp,Tfeed,Tbp)))+vapMixFeed)./vapMixFeed;
qlineSlope = -q/(1-q);


for i=1:length(data)
    t = data(i,3);
    Hmeth = integral(Mcp,298.15,(t+273.15));
    Hiso = integral(Icp,298.15,(t+273.15));
    hmix(i) = x(i)*Hmeth+(1-x(i))*Hiso;
    vapMix(i) = x(i)*HvapM+(1-x(i))*HvapI;
    Hmix(i)=hmix(i)+vapMix(i);
end 

figure;
%hxy
plot(x,hmix./1000,y,Hmix./1000),xlabel("x/y Composition Methanol"),...
    ylabel("Molar Enthalpy (kJ/mol)"),title("Hxy Diagram Methanol/Isopropanol")

results = [D B xM alphaMI slopeROL qlineSlope];
vars = ["Distillate" "Bottoms" "Mole fraction M" "Relative Volatility M/I" "Slope of Rectifying Line" "Slope of Feed Line"];
X=[vars;results];
xlswrite('PrelabSolutions.xlsx',X)
