function [tmpQ] = snowQ_alg(m,a,tsSnow,b,c,tsQshort,d,tsQlong);
% This is  the algorithm that calculates Q at each time step.

% This code is explaiend in detail in: 
% Sproles, E.A., Kerr, T., Orrego Nelson, C. et al. Water Resour Management
% (2016) 30: 2581. doi:10.1007/s11269-016-1271-4

    tmpQ = a*((tsSnow(m-1)).^b) + c*tsQshort(m-1) + d*(tsQlong(m-1));

end