function Dylan_Task1

%{
Name: Dylan Hematillake
ID: 20651646

Program Description: 

Compute the feed temperature and compositions of a given feed for a 
flash drum seperation using the Peng Robinson EOS

%}

clear all
clc
 
global T P R n z
 
R=0.08206; % gas constant.
 
T=50; % Temperature in oC
Pop= [7 10 13 16 19 22 25];
num=length(Pop);
for i=1:num
    P=Pop(i);
    
    F=1; % Feed molar flow rate, mol/s.
    z=[0.3977 0.2926 0.1997 0.0731 0.0369]; % Feed composition in mole fraction.
    n=length(z);
 
    mo=[0.2 0.1 0.2 0.2 0.3 0.5 0.3 0.1 0.05 0.05 0.5]; %first five values are liquid mole fractions, next 5 are vapour, last value is alpha
  
    m=fsolve(@ProcessEquations,mo);
 
    x=m(1:n); % Mole fractions of liquid phase.
	y=m(n+1:2*n); % Mole fraction of vapor phase.
    alpha=m(2*n+1) ; % Ratio of vaporization.
    Results(i,:) = [Pop(i) x y alpha];
end

xlswrite('Dylan Hematillake Lab 6.xlsx',Results,1)

end

function f=ProcessEquations(m)
 
global T P n z
 
x=m(1:n);  % Mole fractions of liquid phase.
y=m(n+1:2*n); % Mole fraction of vapor phase.
alpha=m(2*n+1); % Ratio of vaporization.
 
[phiL,Zm]=fugacity(x,0);
[phiG,Zm]=fugacity(y,1);
 
k=phiL./phiG;
 
f(1:n)=y'-k.*x'; %Equation 5
f(n+1:2*n)= x'-z'./(1+alpha.*(k-1)); %Equation 6(make sure everything is equal to 0                                           
f(2*n+1)= sum(z'.*(1-k)./(1+alpha.*(k-1)));  %Equation 8                                           
 
end

function [phi,Zm]=fugacity(mf,W)
 
global T P R  
 
% Critical values:
CPT=[48.2  305.5  0.099;  % Ethane
     42.0  370.0  0.153;  % Propane
     37.5  425.2  0.199;  % Butane
     33.75 469.60 0.254;  % Pentane
     30.32 507.90 0.300]; % Hexane
 
a=0.45724.*R.^2.*CPT(:,2).^2./CPT(:,1).*(1+(0.37464+1.54226.*CPT(:,3)-0.26992.*CPT(:,3).^2).*(1-sqrt((T+273.15)./CPT(:,2)))).^2;                                                                                      % Equation 12
b=0.0778.*R.*CPT(:,2)./CPT(:,1);  % Equation 13
 
A=a*P/R^2/(T+273.15)^2;
B=b*P/R/(T+273.15);
 
Aij=(1-0.012)*(A*A').^0.5;
Am=mf*(Aij*mf');
Bm=mf*B;
poly=[1 Bm-1 Am-2*Bm-3*Bm^2 Bm^3+Bm^2-Am*Bm]; %Equation 14                                                 
Z=roots(poly);
if W==0
 Zm=min(Z);
else Zm=max(Z);
end
 
phi=exp(B./Bm.*(Zm-1)-log(Zm-Bm)-Am./2/sqrt(2)./Bm.*(2.*(Aij*mf')/Am-B./Bm).*log((Zm+(1+sqrt(2)).*Bm)./(Zm+(1-sqrt(2)).*Bm))); %Equation 27
                                                                                                                                       
end 
 
