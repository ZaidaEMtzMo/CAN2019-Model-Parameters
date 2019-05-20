N = length(obs); x = zeros(N,Ds);
Iinj_ = I_total(1:N)*1e-3;

S0 = (Xh_s(1:Ds,1:N))';
S = (Xh_s(1:Ds,1:N))';
param.gNA = w0_(1); param.gK = w0_(2); param.gL = w0_(3);
% param.gNA = w_(1); param.gK = w_(2); param.gL = w_(3);

obs_ = obs(1:N);
x0 = S0; % the state vector comprising [x1,x2,x3,x4,x5,x6] for All the time samples
x = S;% zeros(size(x0));

dt = param.dt;
A = param.A * 1e-8; Beta_m = param.Beta_m; Gama_m = param.Gama_m; Beta_n = param.Beta_n; Gama_n = param.Gama_n; Phi = param.Phi; C = param.C;
EL = param.EL; ENA = param.ENA; EK = param.EK;

gL = param.gL;  % mS / cm^2 
gNA = param.gNA;    
gK = param.gK;

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
dvdt = (Iinj_*1e-6 / A - I_L - I_K - I_NA) / (C * 1e-3);
% --- *** Very Important to make an accurate estimate of the output
x(2:end,1) = x0(1:end-1,1) + dt*dvdt(2:end);
x(2:end,2) = x0(1:end-1,2) + dt*dndt(2:end);

%         x(:,4) = x0(:,4)+ dt* (A*a*(Ii_        +(vm./(1+exp(r*(v0-x0(:,2)+x0(:,3))))))   -(2*a*x0(:,4))-(a^2*x0(:,1)));
%         x(:,5) = x0(:,5)+ dt* (kA*ka*A*a*(mu_+Ip_+(C2*vm./(1+exp(r*(v0-C1*x0(:,1))))))       -(2*ka*a*x0(:,5))-(ka*ka*a^2*x0(:,2)));
%         x(:,6) = x0(:,6)+ dt* (B*b*(Ip_         +(C4*vm./(1+exp(r*(v0-C3*x0(:,1))))))        -(2*b*x0(:,6))-(b^2*x0(:,3)));        
%         x(:,1) = x0(:,1)+x0(:,4)*dt;
% %         x(:,2) = x0(:,2)+x(:,5)*dt;
% %         x(:,3) = x0(:,3)+x(:,6)*dt;
%         
%         
%         x(1,2) = x0(1,2)+x0(1,5)*dt; 
%         x(1,3) = x0(1,2)+x0(1,6)*dt;
%         
%         x(2:end,2) = x0(2:end,2) + x(1:end-1,5)*dt;% 
%         x(2:end,3) = x0(2:end,3) + x(1:end-1,6)*dt;%


outp = x(:,1);
%outp = outp - mean(outp);
F = norm((obs_- outp),2)

figure; hold on,
plot(yy(:,1),'k')%plot(obs_,'k')
plot(outp,'b')