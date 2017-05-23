function [estSim, param_value, est_nse, val_nse, goodQval, goodQest, goodQval2] = SCR_model(Q, N, nr);
tic

% This Snow Cover Runoff (SCR) model uses monthly calculations of snow cover
% frequency (SCF) from the Modis instrument along with previous monthly streamflow to predict
% monthly streamflow with a one month lead time.

%% The SCR model is based upon the paper:
% Sproles, E.A., Kerr, T., Orrego Nelson, C. Water Resources Management
% (2016) 30: 2581. doi:10.1007/s11269-016-1271-4


% SCF is calculated using Google earth Engine,
% and the code can be found on this same github site. 

% This sample code is set up to run for the La Laguna sub-watershed of the r?o Elqui
% in northern central Chile. The model has also been tested on the John Day
% (easten Oregon, USA) and the r?o Aragon in northern Spain.

%% Model structure

% This code should be cited as:


% The model has seven parameters that are calibrated using this model based
% upon Monte Carlo simulations and Dotty Plots.

% The data that accompanies this code has four columns (Year, Month, Q , N)    

%The inputs are
% Q is a monthly time series of stream flow
% N is the snowcover frequency for the basin, and should be a decimal value.
% nr is the number of realizations for the monte carlo simulation. A good
% number to start with is 1000. 

% The model has seven parameters based upon the formula:
% Qpredicted = a*((tsSnow(m-1)).^b) + c*tsQshort(m-1) + d*(tsQlong(m-1));

% tsSnow is the moving timeseries for SCF. The model optimizes the # of
    % months for the time series.
    
% tsQshort is the moving timeseries for Q in the short term. The model optimizes the # of
    % months for the time series.
    
% tsQlong is the moving timeseries for Q in the long term, and conceptually reresents baseflowconditions
    % The model optimizes the # of months for the time series.
    
% a, c, and d are scaling coefficients specific to each model forcing.
% b is an exponential scaling parameter and represents the tapering effect of snowpack contributions to streamflow as it melts (Leibowitz et al. 2012). 

% this program requires 3 subroutines to be in the same directory:
    % dotty_monte_carsnow.m - creates dotty plots of the calibration realizations.
    % snowQ_alg.m - the actual model algorithm that runs as a subroutine.
    % nashsutcliffe.m - calculates the Nash sutcliffe efficiency, and was
    % developed by Mark Raleigh.
    % additioanlly this program requires MATLAB's Financial Toolbox for the time series
    
%% the actual model...
startDate = datenum('07-01-2002'); % date that time series starts
endDate = datenum('03-01-2017'); % date that time series ends

calWindow_begin = 69; % This starts the steps where the model calibration begins
calWindow_end = 128; %  This starts the steps where the model calibration ends
% The calibration window is 60 months 
% The calibration windows will be different with each basin. The user
% should customize accordingly. 

val = (129:length(N)); % This corresponds with a validation window (48 months)
val2 = (9:68); % This corresponds with a second validation window (60 months)

% Grab the data for the model validation at the end
Qval1 = Q(val);
Qval2 = Q(val2); 
QinNSE = [Qval1;Qval2];

%% Start the Monte Carlo Simulation
% the premise of a Monte Carlo simulation is to randomly select a paramter
% value over a large number of realizations. Then begin to narrow down
% (optimize) the parameters values.

param_value = zeros(7,nr);
est_nse = zeros(1,nr); % nse (Nash-Sutcliffe Efficiency) for the calibration period for the each realization 
val_nse = zeros(1,nr); % nse  for the validation period for the each realization 

h = waitbar(0,'Chugging along...');
for i=1:nr
    waitbar(i/nr)
        % define the lag time...
            % define the values for window_snow (this will become tsSnow)
            hi=12;
            lo=1;
                 window_snow=round(lo+(hi-lo)*rand(1));

            % define the values for window_Qshort (this will become Qshort)
            hi=6;
            lo=2;
                 window_Qshort=round(lo+(hi-lo)*rand(1));
            
            % define the values for window_Qlong (this will become Qlong)
            hi=12;
            lo=7;
                 window_Qlong=round(lo+(hi-lo)*rand(1));
    % Calculate the moving averages. This requires the Financial Toolbox in
    % MATLAB
    tsSnow = tsmovavg(N,'s',window_snow,1);
    tsQshort = tsmovavg(Q,'s',window_Qshort,1);
    tsQlong = tsmovavg(Q,'s',window_Qlong,1);
    
    % Grab the snow data for each realization
    Nest = tsSnow(calWindow_begin:calWindow_end);
    Nval = tsSnow(val); 
    Nval2 = tsSnow(val2); 
    
    % Grab the Qshort data for each realization
    QSest = tsQshort(calWindow_begin:calWindow_end);
    QSval = tsQshort(val);  
    QSval2 = tsQshort(val2);
    
    % Grab the Qlong data for each realization
    QLest = tsQlong(calWindow_begin:calWindow_end);
    QLval = tsQlong(val);
    QLval2 = tsQlong(val2);
    
    % define the parameter values
        % define the value for a
        hi=12;
        lo=0.1;
             a=(lo+(hi-lo)*rand(1));     

        % define the value for b
        hi=20;
        lo=0.1;
             b=(lo+(hi-lo)*rand(1));
              
  
        % define the value for c
        hi=20;
        lo=-20;
              c=(lo+(hi-lo)*rand(1));

        % define the value for d
        hi=10;
        lo=-10;
             d=(lo+(hi-lo)*rand(1));
 
       % Place the parameter values in a matrix to use in the model      
       param_value(1,i) = window_snow;
       param_value(2,i) = window_Qshort;
       param_value(3,i) = window_Qlong;
       param_value(4,i) = a;
       param_value(5,i) = b;
       param_value(6,i) = c;
       param_value(7,i) = d;
        
       holdQest = zeros(1,length(Nest)); % set up a matrix for output values
       
       % Run the model for all of the Monte Carlo Simulations (snowQ_alg is
       % the actual snow cover runoff model)
    
       for m = 2:length(Nest)
           [tmpQest] = snowQ_alg(m,a,Nest,b,c,QSest,d,QLest);
           holdQest(m) = tmpQest;
       end
           holdQest(1) = holdQest(2); % the first month is blank, so fill it with month two. 
        
    % now calculate NSE for the model estimation/calibration data set
     estSim = holdQest';
     estSim(1) = estSim(2); % the first month is blank, so fill it with month two. 
     nse = nashsutcliffe(Q(calWindow_begin:calWindow_end),estSim);
     est_nse(i) = nse;
     
    % Run the model for the validation period 1
    holdQval = zeros(1,length(val));
    q = length(val);
    for z = 2:q
        [tmpQval] = snowQ_alg(z,a,Nval,b,c,QSval,d,QLval);
        holdQval(z) = tmpQval;
    end
        holdQval(1) = holdQval(2); % the first month is blank, so fill it with month two. 
    % Run the model for the validation period 2
    holdQval2 = zeros(1,length(val2));
    q = length(val2);
    for z = 2:q
        [tmpQval2] = snowQ_alg(z,a,Nval2,b,c,QSval2,d,QLval2);
       holdQval2(z) = tmpQval2;
    end
       holdQval2(1) = holdQval2(2); % the first month is blank, so fill it with month two. 
    
    % now calculate NSE for the validation data set
        QsimVal = [holdQval holdQval2];
        nse = nashsutcliffe(QinNSE,QsimVal);
        val_nse(i) = nse;  
end
close(h)
    
est_nse(est_nse<0) = -0.01; %change all of the really bad NSE values to -0.01, mainly for graphing purposes. 
figure('Name','Calibration NSE Values','NumberTitle','off')
plot(est_nse,'o')
ylabel('NSE')
alpha(0.5)
% 
val_nse(val_nse<0) = -0.01; %change all of the really bad NSE values to -0.01, mainly for graphing purposes.
figure('Name','Val NSE Values','NumberTitle','off')
plot(val_nse,'ko')
ylabel('NSE')
alpha(0.5)

% find the best data fit and plot it out...
max_val = max(val_nse)
idx = find(val_nse == max_val);
PB = param_value(:,idx);
% Now reccalculate the good mean time series
    tsSnow = tsmovavg(N,'s',PB(1),1);
    tsQshort = tsmovavg(Q,'s',PB(2),1);
    tsQlong = tsmovavg(Q,'s',PB(3),1);
    
    Nest = tsSnow(calWindow_begin:calWindow_end);
    Nval = tsSnow(val); 
    Nval2 = tsSnow(val2);
    
    QSest = tsQshort(calWindow_begin:calWindow_end);
    QSval = tsQshort(val);  
    QSval2 = tsQshort(val2); 
    
    QLest = tsQlong(calWindow_begin:calWindow_end);
    QLval = tsQlong(val);
    QLval2 = tsQlong(val2);
    
    a = PB(4);
    b = PB(5);
    c = PB(6);
    d = PB(7);
    
    for m = 2:length(Nest) % Run the SCR with best parameters for the calibration data
        [tmpQest] = snowQ_alg(m,PB(4),Nest,PB(5),PB(6),QSest,PB(7),QLest);
        goodQest(m) = tmpQest;
    end
    goodQest(1) = goodQest(2);
    
q = length(val);
    for z = 2:q % Run the SCR with best parameters for the validation data
       [tmpQval] = snowQ_alg(z,PB(4),Nval,PB(5),PB(6),QSval,PB(7),QLval);
       goodQval(z) = tmpQval;
    end
    goodQval(1) = goodQval(2);
    
q = length(val2);
    for z = 2:q
       [tmpQval2] = snowQ_alg(z,PB(4),Nval2,PB(5),PB(6),QSval2,PB(7),QLval2);
       goodQval2(z) = tmpQval2;
    end
    goodQval2(1) = goodQval2(2);
    
    
     dotty_monte_carsnow(est_nse,param_value)
% Plot out the time series
    figure('Name','The realization with the Best Fit','NumberTitle','off','Position', [100, 100, 1049, 895])
    xData = linspace(startDate, endDate, 177);
    plot(xData(69:128), goodQest,'rs')
        hold on
    plot(xData(129:177), goodQval,'bo')
    legend('Calibration Data','Validation Data')
        hold on
    plot(xData(9:68), goodQval2,'bo')
        hold on
    plot(xData,Q, 'k')
    datetick('x','yyyy')
    ylabel('Q (m^3/s)')
      
     
beep
toc

