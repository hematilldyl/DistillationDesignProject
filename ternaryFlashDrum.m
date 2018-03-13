%Flash Drum Script for Ternary Mixture (Part 1, VI) of design project
%Soving acetone, n-butanol, ethanol system, can be parameterized
%function flashDrum = processParameters(Tb1,Tb2,Tb3,A,B,C,MW1,MW2,MW3)

for T = 56:0.01:117.7
%calculate vapor pressures with Antoine's equation
    Pa = (10^(7.11714-(1210.595/(T+229.664))));
    Pb = (10^(7.36366-(1305.195/(T+173.927))));
    Pe = (10^(8.11220-(1592.684/(T+226.184))));
    %relative volatilities
    alpha = Pa/Pb;
    alpha2 = Pa/Pe
    %material balances
    FiVi = 1.25
    FVj = alpha*(FiVi-1)+1;
    Vi= 400/1.25;
    V=Vi+Vj + Vk
    L=1000-V
    Xb=(200-Vj)/L
    Xa =(400-Vi)/L
    Xe = (400-Vk)/L
    Ya = Vi/V
    Yb = Vj/V
    Ye = Vk/V
    %Calculated K
    Ka = Ya/Xa;
    Kb = Yb/Xb;
    Ke = Ye/Xe;
    KR = 1/(Xa*(Ka/Ke)+ Xe + Xb*(Kb/Ke));
    Pvr = KR*760;
    %Calculate estmation temp based on above algorithm
    Tguess = ((-1/((log10(Pvr)-8.11220)))*1592.684)-226.184;
      if Tguess < T+0.5 && Tguess > T-0.5
           T_flash = Tguess;
           break;
      end 
end    
%Material Balance based on design desires   
%Feed basis: 1000 lbmol/h
A = [Ya Xa;Yb Xb; Ye Xe];
b = [400;200;400];
s= A\b;
fprintf('Process parameters: Pressure = 1 atm, T_flash = %.2f %s\n',T_flash,'C')
disp('Feed:  1000 lbmol/h, Acetone: 400 lbmol/h, 1-Butanol: 200 lbmol/h, Ethanol: 400 lbmol/h');
fprintf('Vapour: %.2f lbmol/h, Acetone: %.2f lbmol/h, 1-Butanol: %.2f lbmol/h, Ethanol: %.2f lbmol/h',s(2),Ya*s(2),Yb*s(2), Ye*s(2))
fprintf('\nLiquid: %.2f lbmol/h, Acetone: %.2f lbmol/h, 1-Butanol: %.2f lbmol/h, Ethanol: %.2f lbmol/h',s(1),Xa*s(1),Xb*s(1),Xe*s(1))

%drum sizing
%calculate densities
rhoL = Xa* 0.791 + Xb*0.81 + 0.789*Xe;
rhoV = ((Ya*58.08+Yb*74.121+Ye*46.044)
/(82.0575*(T_flash+273.15)));
%molecular weight of mixture estimation
MWL = (74.121*Xb + 58.08*Xa + 46.044*Xe);
MWV = (Ya*58.08+Yb*74.121 + Ye*46.044);
%material balance
F = ((s(1)*MWL)/(s(2)*MWV))*sqrt(rhoV/rhoL);
%K parameter for the drum (empirical equation)
Kdrum = exp(-1.877478097-0.8145804597*log(F)-0.1870744085*((log(F))^2)-0.0145228667*((log(F))^3)-0.0010148518*((log(F))^4));
%perm velocity in the drum
Uperm = Kdrum*sqrt((rhoL-rhoV)/rhoV);
%cross sectional area
Ac = ((s(2)*MWV*454)/(Uperm*3600*rhoV*28316.85));
Ac = ceil(Ac)+0.5;
D = ceil(sqrt((4*Ac)/pi));
%find height by rule of thumb for vertical flash drum
ht= D*4;


