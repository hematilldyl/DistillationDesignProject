function DylanTask13

%{
Author: Dylan Hematillake
ID: 20651646
Date: 5/16/2018
%}

    clear 
    clc
    format long

    %import data
    time=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'A3:A75');
    Conc=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'B3:B75');

    %inital variables

    CAo=0.050;
    CBo=0.041; %mol/L
    h = 0.5;
    n=length(Conc);

    for i = 1:1:n-4
        CAf1(i) = Conc(i+2);
        tf1(i) = time(i+2);
        CBf1(i) = CBo-(CAo-CAf1(i));
        rf1(i)=-(-3*Conc(i+2)+4*Conc(i+3)-Conc(i+4))/(2*h);
    end

    % Calculate reaction rates based on the forward difference equation using array operations
    global CA CB

    t=time(3:n-2);
    CA=Conc(3:n-2);
    CB=CBo-(CAo-CA);
    rForward=(3*Conc(3:n-2)-4*Conc(4:n-1)+Conc(5:n))./(2*h);

    %central difference formula
    rCentral = -(Conc(4:n-1)-(Conc(2:n-3)))/(2*h);

    %gradient
    rg = -gradient(CA,t);
    ng = length(rg);

    subplot(1,3,1);

    plot(tf1,rf1,'d',t,rForward,'o',t,rCentral,'*',ng,rg,'r'),xlabel('t'),ylabel('Rate')
    legend('for loop','Forward Difference','Central Difference','Gradient')
    hold on

    %task 3

    m0 = [6 1 1]; %iv
    lb = [0.1 0.01 0.01]; %lower bound
    ub = [10 2 2] ;%upper bound

    options =optimset('TolX',1e-16,'TolFun',1e-16,'MaxFunEval',4000,'MaxIter',4000);
    r = rCentral;
    [m,residn, resid,exiflage,output,lamda,J]=lsqcurvefit(@KineticEquation,m0,t,r,lb,ub,options);
    m=m

    ci = nlparci(m,resid,J) % 95% confidence interval

    subplot(1,3,2)
    plot(t,resid,'o'),xlabel('t'),ylabel("Residual");


    rpred = KineticEquation(m,t); %subfun
    
    subplot(1,3,3)
    plot(t,rCentral,'o',t,rpred,'-'),xlabel('t'),ylabel('r')
    figure;
    plot(t,r-rpred,'o'),xlabel('t'),ylabel('Rate of Reaction')
    SSE = sum((r-rpred).^2);
    rmean = sum(r)/length(r);
    SST = sum((r-rmean).^2);
    R2 = 1-SSE/SST;
    
end

function rpred=KineticEquation(m,t)

    global CA CB

    k = m(1);
    alpha = m(2);
    beta = m(3);

    rpred = k.*CA.^alpha.*CB.^beta;


end 








