function GroupOTask5

    %{
    Authors: Dylan Hematillake, Aaron Bleasdale-Pollowy, 
             Matt Williams, Victoria Bodnitski
    ID: 20651646,,,20657294,
    Date: 5/22/2018
    %}

    clear
    clc

    global CAo CBo CA CB time Conc

    %load data
    time=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'A3:A75');
    Conc=xlsread('Experiment5_KineticDataFromChE101.xlsx', 'B3:B75');

    %starting parameters
    n = length(time); 
    CAo=0.050;
    CBo = 0.041;
    CA = Conc(3:n-2);
    CB = CBo-(CAo-CA);
    h = 0.5;
    N = length(time);


    options = optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEval',2000,'MaxIter',2000);

    m0 = [6 1 1]; %iv
    lb = [0.1 0.01 0.01]; %lower bound
    ub = [10 2 2] ;%upper bound

    %regression
    [m,rnorm,resid,exitflag,output,lamda,J]=lsqcurvefit(@Modeldata,m0,time,Conc,lb,ub,options);
    m=m
    %calculate CI
    ci=nlparci(m,resid,J)

    %plot residuals
    subplot(1,2,1);
    plot(time,resid,'o'),xlabel('t'),ylabel('Residual');

    %calculate predicted values
    CApred=Modeldata(m,time);

    %calculate statistics
    SSE = sum((Conc-CApred).^2)
    Cmean = sum(Conc)/length(Conc)
    SST = sum((Conc-Cmean).^2)
    R2 = 1-SSE/SST

    %plot actual & predicted
    subplot(1,2,2);
    plot(time,Conc,'d',time,CApred,'-'),xlabel('t'),ylabel('r');
    legend("Actual","Predicted")
    
    res = Conc - CApred;
    figure;
    plot(time,res,'o'),xlabel("t"),ylabel("Concentration A")

end

%function model
function CAm = Modeldata(m,t)   
    global CAo CBo time
    
    %initial values
    y0 = [CAo;CBo];
    [t c] = ode45(@ODEEquations,time,y0,[],m);
    CAm = c(:,1);
end


%system of equations for model
function f = ODEEquations(t,C,m)
    
    k = m(1);
    alpha = m(2);
    beta = m(3);
    
    f(1) = -k*(C(1)^alpha)*(C(2)^beta);
    f(2) = -k*(C(1)^alpha)*(C(2)^beta);

    f = f';
    
end 
