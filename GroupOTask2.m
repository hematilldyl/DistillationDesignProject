function GroupOTask2

    %{
    Authors: Dylan Hematillake, Aaron Bleasdale-Pollowy, 
             Matt Williams, Victoria Bodnitski
    ID: 20651646,,,20657294,
    Date: 5/22/2018
    %}

    clear
    clc

    global CA CB

    %load data
    time=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'A3:A75');
    Conc=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'B3:B75');

    %starting paramaters
    n = length(time); 
    CAo=0.050;
    CBo = 0.041;
    CA = Conc(3:n-2);
    CB = CBo-(CAo-CA);
    h = 0.5;
    rc = (Conc(2:n-3)-Conc(4:n-1))/(2*h);
    t=time(3:n-2);
    %linearize equation
    y(:,1) = log(rc);
    X(:,1) = log(CA);
    X(:,2) = log(CB);

    %fit data with regression
    md = fitlm(X,y);
    b = md.Coefficients.Estimate
    bCI = coefCI(md);

    %extract model parameters from transformation
    k = exp(b(1))
    kint = exp(bCI(1,:));
    alpha = b(2)
    beta = b(3)

    %evaluate regression for each point
    ypred = feval(md,X);
    res = y - ypred
    plot(t,res,'o'),xlabel('t'),ylabel('Rate of Reaction')
    %calculate the statistics
    SSE = sum((y-ypred).^2)
    ymean = sum(y)/length(y)
    SST = sum((y-ymean).^2)
    R2 = 1-SSE/SST
    
end
