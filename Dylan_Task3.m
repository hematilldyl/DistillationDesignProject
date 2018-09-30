function Dylan_Task3

%Dylan Hematillake
%20651646

clear all
clc
 
global T P R n z
 
R=0.08206; % gas constant.
 

Pop= [7 10 13 16 19 22 25];
num=length(Pop);
for i=1:num
    P=Pop(i);
    
    F=1; % Feed molar flow rate, mol/s.
    z=[0.3977 0.2926 0.1997 0.0731 0.0369]; % Feed composition in mole fraction.
    n=length(z);
    T0=70 
    mo=[0.1 0.4 0.4 0.4 0.4 T0]; 
  
    m=fsolve(@ProcessEquations,mo);
 
    y=z; %y=z because dew point calculation
    x=m(1:n) ;
    T=m(n+1)  ;
    Results(i,:) = [Pop(i) T];
end

xlswrite('Dylan Hematillake Lab 6.xlsx',Results,3)

end

function f=ProcessEquations(m)
 
global P n z
 
 
y=z; %y=z because dew point calculation
x=m(1:n); 
T=m(n+1); 
 
[phiL,Zm]=fugacity(x,T,0);
[phiG,Zm]=fugacity(y,T,1);
 
k=phiL./phiG;
 
f(1:n)=y'-k.*x'; %Equation 5 
f(n+1)= sum(z'./k)-1; %Equation 10                                        
 
end

function [phi,Zm]=fugacity(mf,T,W)
 
global P R  
 
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
 
