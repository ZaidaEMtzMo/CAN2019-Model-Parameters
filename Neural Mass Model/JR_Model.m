function F = JR_Model(S,mu,Ip,Ii,p)

x0 = S; % the state vector comprising [x1,x2,x3,x4,x5,x6]
x = zeros(size(x0));
dt = p.dt;
A = p.A; a = p.a;
B = p.B; b = p.b;
vm = p.vm; r = p.r;
v0 = p.v0;
C1 = p.C1; C2 = p.C2; C3 = p.C3; C4 = p.C4;
ka = p.ka; kA = p.kA;

        x(1) = x0(1)+x0(4)*dt;
        x(2) = x0(2)+x0(5)*dt;
        x(3) = x0(3)+x0(6)*dt;
        x(4) = x0(4)+ dt* (A*a*(Ii        +(vm/(1+exp(r*(v0-x0(2)+x0(3))))))   -(2*a*x0(4))-(a^2*x0(1)));
        x(5) = x0(5)+ dt* (kA*ka*A*a*(mu(1)+Ip+(C2*vm/(1+exp(r*(v0-C1*x0(1))))))       -(2*ka*a*x0(5))-(ka*ka*a^2*x0(2)));
        x(6) = x0(6)+ dt* (B*b*(Ip         +(C4*vm/(1+exp(r*(v0-C3*x0(1))))))        -(2*b*x0(6))-(b^2*x0(3)));

F = x;

