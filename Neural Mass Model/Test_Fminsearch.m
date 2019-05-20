N = length(indx);
obs_ = obs(indx(2):indx(N));
S0 = (Xh_s(1:6,indx(1):indx(N-1)))';
S = (Xh_s(1:6,indx(2):indx(N)))';
fun = @(x)Optimization_JR_Model(obs_,S0,S,x,mu(indx(1:end-1)),Ip(indx(1:end-1))',Ii(indx(1:end-1))',p);

% x0_ = [500,500,500,500];
options = optimset('PlotFcns',@optimplotfval);
% run test_FullUnderstanding_ParameterEstimation
[x_,fval] = fminsearch(fun,x0_,options);