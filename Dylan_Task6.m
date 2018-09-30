function Dylan_Task6

%Dylan Hematillake
%20651646

clear all
clc

global P x z
Tguess = [133 120 110]
Pop = [300 100 50]./14.69594878
Title = ["Txy 300 psia","Txy 100 psia","Txy 50 psia","xy 300 psia", "xy 100 psia".....
        ,"xy 50 psia"];
index = 4;
for j = 1:3
    for i = 1:101
        P = Pop(j);
        xb(i)=(i-1)*0.01;
        x=[xb(i) 1-xb(i)];
        z = [0.5 0.5]';

        yb0 = (i-1)*0.01;
        T0 = Tguess(j)-(110-33)/101*i;
        alpha0 = (i-1)*0.01;
        mo = [yb0 T0 alpha0];

        m = fsolve(@ProcessEquations, mo);
        yb(i)=m(1);
        Tb(i) = m(2);
        alpha(i)=m(3);
        
        
    end

    Tb = 273.15+Tb;
    M = [Tb',xb',yb'];
    figure;
    subplot(1,2,1);
    plot(xb,Tb,yb,Tb);title(Title(j)),xlabel("x/y"),ylabel("Temperature (K)"),......
        xlim([0,1])
    subplot(1,2,2);
    plot(xb,yb,xb,xb);title(Title(index)),xlabel("x"),ylabel("y")
    index=index+1;
    xlswrite('Task6.xlsx',M,j)
end

end

function f = ProcessEquations(m)

global x z 

y = [m(1) 1-m(1)];
T=m(2);
alpha=m(3);

phiL = fugacity(x,T,0);
phiG = fugacity(y,T,1);

k = phiL./phiG;

f(1:2) = y'-k.*x';
f(3)=sum(z.*(1-k)./(1+alpha.*(k-1)));

end
function [phi,Zm]=fugacity(mf,T,W)
 
global P R  

R = 0.08206;
% Critical values:
CPT=[42.0  370.0  0.153;  % Propane
     37.5  425.2  0.199;]; %butane
 

a=0.45724.*R.^2.*CPT(:,2).^2./CPT(:,1).*(1+(0.37464+1.54226.*CPT(:,3)-0.26992.*CPT(:,3).^2).*(1-sqrt((T+273.15)./CPT(:,2)))).^2;                                                                                      % Equation 12
b=0.0778.*R.*CPT(:,2)./CPT(:,1);  % Equation 13
 
A=a*P/R^2/(T+273.15)^2;
B=b*P/R/(T+273.15);
 
Aij=(A*A').^0.5;
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



    
    
