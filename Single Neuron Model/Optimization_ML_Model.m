function F = Optimization_ML_Model(obs,S,init_pp,I_inj_,pp)

N = length(obs); x = zeros(N,2);
x0 = S;
gNA = init_pp(1); gK = init_pp(2); gL = init_pp(3);

dt = pp.dt;
A = pp.A * 1e-8; Beta_m = pp.Beta_m; Gama_m = pp.Gama_m; Beta_n = pp.Beta_n; Gama_n = pp.Gama_n; Phi = pp.Phi; C = pp.C;
EL = pp.EL; ENA = pp.ENA; EK = pp.EK;

%--- Leak current
I_L = gL * (x0(:,1) - EL) * 1e-3; % mA / cm^2
%--- Sodium current and m_inf
minf = 0.5*(1+tanh((x0(:,1) - Beta_m)/Gama_m));
I_NA = gNA * minf.*(x0(:,1) - ENA) * 1e-3;  % mA / cm^2
%--- Potassium current
nInf = 0.5*(1+tanh((x0(:,1) - Beta_n)/Gama_n));%alpha_n(j)/(alpha_n(j) + beta_n(j));
taun = 1./cosh(0.5*(x0(:,1) - Beta_n)/Beta_n);
I_K = gK*x0(:,2).*(x0(:,1) - EK) * 1e-3;  % mA / cm^2  

%--- membrane potential
dndt = (nInf - x0(:,2))*Phi ./ taun;
dvdt = (I_inj_*1e-6 / A - I_L - I_K - I_NA) / (C * 1e-3);
% --- *** Very Important to make an accurate estimate of the output
x(2:end,1) = x0(1:end-1,1) + dt*dvdt(2:end);
x(2:end,2) = x0(1:end-1,2) + dt*dndt(2:end);

outp = x(:,1);
%outp = outp - mean(outp);
F = norm((obs - outp),2);