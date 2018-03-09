function ret= txypxy(A,B,C,P,Tp, Species1,Species2)
%{
Author: Dylan Hematillake
Date: April 7th, 2017
Program Description:
    Function takes Antoine's Constants for T in Celcius and P in mmHg, A B C in vector form with Species 1 at
    index 1, Species 2 at index 2. Then the constant pressure for the txy
    diagram is taken as well as the constant temperature for the pxy diagram.
    The names of species 1 and 2 in string format are recieved as well and the
    plots are produced.
%}

%Variable declarations
Txypressure = num2str(P);
temperaturePxy = num2str(Tp);
Tb2 = (B(2)/(A(2)-log10(P))) - C(2);
Tb1 = (B(1)/(A(1)-log10(P))) - C(1);
x1(1) = 0;

T = Tb2:-1:Tb1;
for i = 1:length(T)
    P1(i) = 10^(A(1)-(B(1)/(C(1)+T(i))));
    P2(i) = 10^(A(2)-(B(2)/(C(2)+T(i))));
    if i ~=1
    x1(i) = (P-P2(i))/(P1(i)-P2(i));
    end
    y1(i)= (x1(i)*P1(i))/P;
end
x1(length(x1))=1;
y1(length(y1))=1;
%txy generation
figure;
subplot(2,1,1);
plot(x1,T,y1,T);
title(['Txy Diagram for ' Species1 ' and ' Species2 ' at ' Txypressure ' mmHg']);
xlabel(['Mole fraction of ' Species1]);
ylabel('Temperature (C)');

%pxy generation 
xi = 0:0.01:1;
Pv1= 10^(A(1)-(B(1)/(C(1)+Tp)));
Pv2 = 10^(A(2)-(B(2)/(C(2)+Tp)));
Pt(1) = Pv2;
for j = 1:100
    if (j ~= 1 || j~=100) 
    Pt(j) = xi(j)*Pv1+(1-xi(j))*Pv2;
    end
    yi(j) = (xi(j)*Pv1)/Pt(j);
end
Pt(length(Pt)+1) = Pv1;
yi(length(yi)+1) = (xi(length(xi))*Pv1)/Pt(length(Pt));

subplot(2,1,2);
plot(xi,Pt,yi,Pt);
title(['Pxy Diagram for ' Species1 ' and ' Species2 ' at ' temperaturePxy ' Celcius']);
xlabel(['Mole fraction of ' Species1]);
ylabel('Pressure (mmHg)');
end
