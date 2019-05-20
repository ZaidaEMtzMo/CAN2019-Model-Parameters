function F = Optimization_JR_Model(obs,S0,S,init_param,mu,Ip,Ii,p)

C1 = init_param(1); C2 = init_param(2); C3 = init_param(3); C4 = init_param(4);

x0 = S0; % the state vector comprising [x1,x2,x3,x4,x5,x6] for All the time samples
x = S;% zeros(size(x0));
dt = p.dt;
A = p.A; a = p.a;
B = p.B; b = p.b;
vm = p.vm; r = p.r;
v0 = p.v0;

ka = p.ka; kA = p.kA;

        x(:,4) = x0(:,4)+ dt* (A*a*(Ii        +(vm./(1+exp(r*(v0-x0(:,2)+x0(:,3))))))   -(2*a*x0(:,4))-(a^2*x0(:,1)));
        x(:,5) = x0(:,5)+ dt* (kA*ka*A*a*(mu+Ip+(C2*vm./(1+exp(r*(v0-C1*x0(:,1))))))       -(2*ka*a*x0(:,5))-(ka*ka*a^2*x0(:,2)));
        x(:,6) = x0(:,6)+ dt* (B*b*(Ip         +(C4*vm./(1+exp(r*(v0-C3*x0(:,1))))))        -(2*b*x0(:,6))-(b^2*x0(:,3)));        
        x(:,1) = x0(:,1)+x0(:,4)*dt;
%         x(:,2) = x0(:,2)+x(:,5)*dt;
%         x(:,3) = x0(:,3)+x(:,6)*dt;

        x(1,2) = x0(1,2)+x0(1,5)*dt; 
        x(1,3) = x0(1,2)+x0(1,6)*dt;
        
        x(2:end,2) = x0(2:end,2) + x(1:end-1,5)*dt; 
        x(2:end,3) = x0(2:end,3) + x(1:end-1,6)*dt;


outp = x(2:end,2) - x(2:end,3);
outp = outp - mean(outp);
obs_ = obs - mean(obs);
F = norm((obs_(2:end) - outp),2);

