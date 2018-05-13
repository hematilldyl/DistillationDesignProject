%Solve for activity coefficients, sat. vapor pressure, vapor compositions and Temperature of binary ideal vapor, non ideal liquid system
%Uses van Laar equations to relate activity coefficients to compositions

function F = roots(x)
%System of equations to solve
global n
 
F(1) = log(x(1))-(2.96./(1+((2.96.*n)/(1.63.*(1-n)))).^2);
F(2) = log(x(2))-(1.63./(1+((1.63.*(1-n))/(2.96.*(n)))).^2);
F(3) = n.*x(1).*x(3)-x(5).*0.709275;
F(4) = (1-n).*x(2).*x(4)-x(6).*0.709275;
F(5) = log10(x(3))-(4.13983-(1316.554./(x(7)-35.581)));
F(6) = log10(x(4))-(5.24677-(1598.673./(x(7)-46.424)));
F(7) = 1-x(5)-x(6);
 
end
_____________________________________________________________________________
clear
clc 
 
 
global n
 
%variables: gC,gE,Pvc,Pve,yc,ye,T
 
%initial guess
x = [1,1,0.1,0.1,0.1,0.1,273.15];
%convert to functions
 
 
X =zeros(9,7);
i = 1;

for n=0.1:0.1:0.9
   
    F = @roots;
    x = fsolve(F,x);
    X(i,:) = x;
    i = i+1;
end

xlswrite('simulated_txy_data.xlsx',X)
