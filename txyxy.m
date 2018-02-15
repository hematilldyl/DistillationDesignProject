%X-y and T-x-y script
function ret = txy(A,B,C,P,Species1,Species2,comp)
%Acetone Parameters: A: 7.11714 B: 1210.595 C: 229.664 Tb: 56
%n-Butanol Parameters: A: 7.36366 B: 1305.198 C: 173.427  Tb=117.40
Txypressure = num2str(P);
Tb2 = (B(2)/(A(2)-log10(P))) - C(2)
Tb1 = (B(1)/(A(1)-log10(P))) - C(1)
if Tb2 >Tb1
    T = Tb2:-0.001:Tb1;
else
    T = Tb1:-0.001:Tb2;
end
for i = 1:length(T)
    P1(i) = 10^(A(1)-(B(1)/(C(1)+T(i))));
    P2(i) = 10^(A(2)-(B(2)/(C(2)+T(i))));
    x1(i) = (P-P2(i))/(P1(i)-P2(i));
    y1(i)= (x1(i)*P1(i))/P;
end
%txy generation
figure;
subplot(2,1,1);
plot(x1,T,y1,T);
title('Txy Diagram for ' +Species1 +' and ' +Species2 +' at ' +Txypressure +' mmHg');
xlabel(['Mole fraction of ' Species1]);
ylabel('Temperature (C)');
%only for acetone and 1-butanol, assume alpha is constant
alpha = exp((41800/8.314)*((1/329.15)-(1/390.85)));
for i = 1:100
    frac = (i-1)*0.01;
    x(i)= frac;
    y(i) = (alpha*frac)/(1+(alpha-1)*frac);
end
subplot(2,1,2); plot(x,y,y,y) ;title('xy Diagram for ' +Species1 +' and ' +Species2 );
xlabel('Liquid mole fraction x'); ylabel('Vapour mole fraction y');
