function DylanTask4

%{
Author: Dylan Hematillake
ID: 20651646
Date: 5/16/2018
%}
    clear all
    clc
    format long

    t=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'A3:A75');
    Conc=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'B3:B75');
    
    CAo = Conc(1);
    CBo = 0.041;
    N = length(t);
    
    t1 = t(1:N-10);
    CA1 = Conc(1:N-10);
    
    y = log(CA1./(CBo-CAo+CA1));
    md = fitlm(t1,y)
    b = md.Coefficients.Estimate
    bCI = coefCI(md)
    kp = b(2)/(CAo-CBo)
    kpCI = bCI(2,:)/(CAo-CBo)
    ypred = feval(md,t1);
    
    figure;
    plot(t1,y,'o',t1,ypred,'-'),xlabel('t'),ylabel('LHS Linearized Rate Equation');
    legend('Experimental Values','Predicted Values');
    res = y-ypred
    figure;
    plot(t1,res,'o'),xlabel('t'),ylabel('LHS Linearized Rate Equation');
    SSE = sum((y-ypred).^2)
    ymean = sum(y)/length(y)
    SST = sum((y-ymean).^2)
    R2 = 1-SSE/SST
    
end
    
