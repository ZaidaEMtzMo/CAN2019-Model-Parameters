function F = ML_Model(S,Iin_nA,pp)


x0 = S; % the state vector only (not the constant parameters) comprising [V,n]
x = zeros(size(x0));
dt = pp.dt;
A = pp.A * 1e-8; Beta_m = pp.Beta_m; Gama_m = pp.Gama_m; Beta_n = pp.Beta_n; Gama_n = pp.Gama_n; Phi = pp.Phi; C = pp.C;
EL = pp.EL; ENA = pp.ENA; EK = pp.EK;

gL = pp.gL;  % mS / cm^2 
gNA = pp.gNA;    
gK = pp.gK;

v = x0(1);
n = x0(2);

%--- Leak current
I_L = gL * (v - EL) * 1e-3; % mA / cm^2
%--- Sodium current and m_inf
minf = 0.5*(1+tanh((v - Beta_m)/Gama_m));
I_NA = gNA * minf*(v - ENA) * 1e-3;  % mA / cm^2
%--- Potassium current
nInf = 0.5*(1+tanh((v - Beta_n)/Gama_n));%alpha_n(j)/(alpha_n(j) + beta_n(j));
taun = 1./cosh(0.5*(v - Beta_n)/Beta_n);
I_K = gK*n*(v - EK) * 1e-3;  % mA / cm^2  

%--- membrane potential
dndt = (nInf - n)*Phi / taun;
dvdt = (Iin_nA*1e-6 / A - I_L - I_K - I_NA) / (C * 1e-3);

x(1) = x0(1) + dt*dvdt;
x(2) = x0(2) + dt*dndt;

F = x;