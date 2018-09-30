function GroupOTask6

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

    %inital parameters
    n = length(time); 
    CAo=0.050;
    CBo = 0.041;
    CA = Conc(3:n-2);
    CB = CBo-(CAo-CA);
    h = 0.5;
    N = length(time);


    options =optimset('TolX',1e-16,'TolFun',1e-16,'MaxFunEval',2000,'MaxIter',2000);

    m0 = [6 100 100]; %iv
    lb = [0.001 10 10]; %lower bound
    ub = [10 500 500] ;%upper bound

    %regression on data
    [m,rnorm,resid,exitflag,output,lamda,J]=lsqcurvefit(@Modeldata,m0,time,Conc,lb,ub,options);
    m=m;
    %calculate CI
    ci=nlparci(m,resid,J)

    %plot residuals
    subplot(1,2,1);
    plot(time,resid,'o'),xlabel('t'),ylabel('Residual');

    %calculate predicted values
    CApred=Modeldata(m,time);
    res = Conc-CApred
    %calculate statistics
    SSE = sum((Conc-CApred).^2);
    Cmean = sum(Conc)/length(Conc);
    SST = sum((Conc-Cmean).^2);
    R2 = 1-SSE/SST;

    %plot actual and predicted values
    subplot(1,2,2);
    plot(time,Conc,'d',time,CApred,'-'),xlabel('t'),ylabel('Concentration');
    legend("Actual","Predicted")
end

%model equation to solve
function CAm = Modeldata(m,t)   
    global CAo CBo time
    
    %initial values
    y0 = [CAo;CBo];
    [t c] = ode45(@ODEEquations,time,y0,[],m);
    CAm = c(:,1);
end

%system of equations for hydrolysis
function f = ODEEquations(t,C,m)
    
    k1 = m(1);
    k2 = m(2);
    k3 = m(3);
    
    f(1) = -k1*(1-1/(k2+k3))*C(1)*C(2);
    f(2) = -k1*(1-1/(k2+k3))*C(1)*C(2);
    f = f';
    
end 
