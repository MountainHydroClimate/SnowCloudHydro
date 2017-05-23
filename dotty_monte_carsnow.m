function dotty_monte_carsnow(realization_nse,param_value)
%This subroutine is called by the Snow Cover Runoff (SCR) model...
% and plots out a series of dotty plots. Each dot corresponds to the parameter value (x-axis) and
% the Nash Sutcliffe Efficiency (realization_nse) for each realization.

figure('Name','NSE Values by Parameter','NumberTitle','off')
param_str = {'Snow Window', 'Q Short Time Window','Q Long Time Window','a','b','c','d'};

for oo = 1:7
    minp(oo) = min(param_value(oo,:)) -0.5;
    maxp(oo) = max(param_value(oo,:)) -0.5;
end

for p = 1:7
    subplot(4,2,p)
    plot(param_value(p,:),realization_nse(:),'k.')
    alpha(0.5)
    hold on
    ylim([0 1]);
    ylabel('NSE')
    xlabel('Parameter Value')
    str_param = (param_str(p));
    title(str_param)
end

