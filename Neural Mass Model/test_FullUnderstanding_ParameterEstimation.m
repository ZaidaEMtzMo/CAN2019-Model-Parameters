N = length(indx); x = zeros(N-1,6);
mu_ = mu(indx(1):indx(N-1)); Ip_ = Ip(indx(1):indx(N-1))'; Ii_ = Ii(indx(1):indx(N-1))';

% C=135; C1=C;
% C2=0.8*C; C3=0.25*C;
% C4=0.25*C;
% C1 = 50; C2 = 60; C3 = 60; C4 = 50;
% C1 = x0_(1); C2 = x0_(2); C3 = x0_(3); C4 = x0_(4);
C1 = x_(1); C2 = x_(2); C3 = x_(3); C4 = x_(4);

S0 = (Xh_s(1:6,indx(1):indx(N-1)))';
S = (Xh_s(1:6,indx(2):indx(N)))';

obs_ = obs(indx(2):indx(N));

x0 = S0; % the state vector comprising [x1,x2,x3,x4,x5,x6] for All the time samples
x = S;% zeros(size(x0));
dt = p.dt;
A = p.A; a = p.a;
B = p.B; b = p.b;
vm = p.vm; r = p.r;
v0 = p.v0;

ka = p.ka; kA = p.kA;

        x(:,4) = x0(:,4)+ dt* (A*a*(Ii_        +(vm./(1+exp(r*(v0-x0(:,2)+x0(:,3))))))   -(2*a*x0(:,4))-(a^2*x0(:,1)));
        x(:,5) = x0(:,5)+ dt* (kA*ka*A*a*(mu_+Ip_+(C2*vm./(1+exp(r*(v0-C1*x0(:,1))))))       -(2*ka*a*x0(:,5))-(ka*ka*a^2*x0(:,2)));
        x(:,6) = x0(:,6)+ dt* (B*b*(Ip_         +(C4*vm./(1+exp(r*(v0-C3*x0(:,1))))))        -(2*b*x0(:,6))-(b^2*x0(:,3)));        
        x(:,1) = x0(:,1)+x0(:,4)*dt;
%         x(:,2) = x0(:,2)+x(:,5)*dt;
%         x(:,3) = x0(:,3)+x(:,6)*dt;
        
        
        x(1,2) = x0(1,2)+x0(1,5)*dt; 
        x(1,3) = x0(1,2)+x0(1,6)*dt;
        
        x(2:end,2) = x0(2:end,2) + x(1:end-1,5)*dt;% 
        x(2:end,3) = x0(2:end,3) + x(1:end-1,6)*dt;%


outp = x(2:end,2) - x(2:end,3);
outp = outp - mean(outp);
F = norm((obs_(1:end-1) - outp),2)

figure; hold on,
plot(obs_(2:end) - mean(obs_(2:end)),'k')
plot(outp,'b')