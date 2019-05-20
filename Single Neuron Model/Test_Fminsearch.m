N = length(tt);
obs_ = obs;
S = (Xh_s(1:Ds,:))';
w0_ = [140,100,20];%[40;50;5];
fun = @(x)Optimization_ML_Model(obs_,S,x,I_total*1e-3,param);
% t_init = 0.5; % sec
% t_init_index = t_init * 1e3/dt;
% w0_ = mean(Xh_s(Ds+1:Ds+Dp,t_init_index:end),2);
% x0_ = [500,500,500,500];
options = optimset('PlotFcns',@optimplotfval);
% run test_FullUnderstanding_ParameterEstimation
[x_,fval] = fminsearch(fun,w0_,options);