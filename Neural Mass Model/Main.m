clear all
%% Parameter Settings
dt=0.001; % sampling time (sec)
F=[10,20,40,80,100,130,160]; % stimulation frequency
stim_time=10000;   % stimulation train duration
Tbase=5000;      % Samples before stimulation starts.
Ttot=16000;       % Total number of samples simulated
% the first row saves 1 node model simulation, 2nd and 3rd row saves 2 node model simulation
sim_lfp=cell(3,length(F));
% Model parameters 
v0=6; vm=5; r=0.3;
C=135; C1=C;
C2=0.8*C; C3=0.25*C;
C4=0.25*C;

A=3.25;
B=22;
a=100;
b=50;
ka=1;kA=1;
%--- define parameter p
p.dt = dt;
p.A = A; p.a = a;
p.B = B; p.b = b;
p.vm = vm; p.r = r; p.v0 = v0;
p.C1 = C1; p.C2 = C2; p.C3 = C3; p.C4 = C4;
p.ka = ka; p.kA = kA;
%% stimulation and running simulation
stim_freq=F(3);
cycle=[1,-1,zeros(1,floor(1000/stim_freq)-2)];
stim=[];
for ss=1:(floor(stim_time*stim_freq/1000))
    stim=[stim,cycle];
end
I=[zeros(1,Tbase),stim,zeros(1,Ttot-Tbase-length(stim))];
Ip=60*I;
Ii=60*r*I;
L=length(I);
X=zeros(L,6);
mu = 0.1*randn(L,1);

for k=2:L
    F = JR_Model(X(k-1,:),mu(k-1),Ip(k-1),Ii(k-1),p);
    X(k,:) = F;
end

out_EEG = X(:,2)-X(:,3);%X(:,2) - X(:,3);% 
out_EEG = out_EEG + 0.01*randn(size(out_EEG));%- mean(out_EEG);% 
tt = 0:dt:(length(out_EEG)-1)*dt;
figure; plot(tt,I)
xlabel('Time (Sec)')
ylabel('Stimulus (normalized)')
title('Stimulus (DBS)')

figure; plot(tt,out_EEG)
xlabel('Time (Sec)')
ylabel('Voltage (mV)')
title('Recorded EEG')
