

Xh_s=[];
Xh_=[];

obs = out_EEG;
L2=1;%size(observation);
Obs_dim=L2;
s_dim = 10; % x1:x6 & C1:C6
nsp_s = 2*s_dim + 1;

%%% Weighths
alpha = 1;
ktau = 0;%
beta = 3;
lambda_s = alpha^2 *(s_dim+ktau) - s_dim;
W_s = [lambda_s 0.5 0]/(s_dim+lambda_s); 
W_s(3) = W_s(1) + (1-alpha^2) + beta;
% W_s(1) = +0;%W_s(1) + (1-alpha^2) + beta;
gamma_s = sqrt(s_dim+lambda_s);
%%%%%%%%%%
Rv_s = 1e-2*eye(s_dim);
Rww_s = 1e-6*eye(Obs_dim);
Sv_s = chol(Rv_s)';
%Sw = chol(Rww_s)';
x1=0; x2=0; x3=0; x4=0; x5=0; x6=0; c1 = 0; c2 = 0; c3 = 0; c4 = 0;
s = [x1;x2;x3;x4;x5;x6;c1;c2;c3;c4];
Xh_s(:,1) = s;
Ps = 1*eye(size(s,1));

aa=0.9;
for j=1:1
for k=1000:length(obs)
 
%Sw = chol(Rww)';   
%%%    State Estimation
    Z=[];
    Sx_s = chol(Ps)';
    Z    = cvecrep(s, nsp_s);                                                   
    sSz  = gamma_s * Sx_s;
    sSzs = [sSz -sSz];
    Z(:,2:nsp_s) = Z(:,2:nsp_s) + sSzs;           % build sigma-point set
    X_=[];
    Y_=[];
    temp1=[];
    for i=1:nsp_s
        p.C1 = Z(7,i); p.C2 = Z(8,i); p.C3 = Z(9,i); p.C4 = Z(10,i);
        aaa = JR_Model(Z(1:6,i),mu(k-1),Ip(k-1),Ii(k-1),p); % F(X_sigma(:,i),input);
        X_(:,i)=[aaa Z(7:10,i)'] + 1e-10*randn(1,10);
    end
    sh_(:,k) = W_s(1)*X_(:,1) + W_s(2)*sum(X_(:,2:nsp_s),2);% + 0.01*randn(size(xh_(:,k-1)));
    temp1 = X_ - cvecrep(sh_(:,k),nsp_s);
    Ps_ = W_s(3)*temp1(:,1)*temp1(:,1)' + W_s(2)*temp1(:,2:nsp_s)*temp1(:,2:nsp_s)'+ Rv_s;
    %%% augmented sigma points
%     Z=[];
%     Sx_s = chol(Ps_)';
%     Z    = cvecrep(sh_(:,k), nsp_s);                                                   
%     sSz  = gamma_s * Sx_s;
%     sSzs = [sSz -sSz];
%     Z(:,2:nsp_s) = Z(:,2:nsp_s) + sSzs;
    %%%
    Y_=[];
    yh_=[];
    temp2=[];
    Y_ = X_(2,:) - X_(3,:);% Y_ = Y_ - mean(Y_); % observation function
    yh_(:,k) = W_s(1)*Y_(:,1) + W_s(2)*sum(Y_(:,2:nsp_s),2);
    temp2 = Y_(1,:) - cvecrep(yh_(1,k),nsp_s);
    Py  = W_s(3)*temp2(:,1)*temp2(:,1)' + W_s(2)*temp2(:,2:nsp_s)*temp2(:,2:nsp_s)' + Rww_s;
   
    Psy = W_s(3)*temp1(:,1)*temp2(:,1)' + W_s(2)*temp1(:,2:nsp_s)*temp2(:,2:nsp_s)';
    KG_s = Psy / Py;
    inov(:,k) = obs(k) - yh_(1,k);
    Xh_s(:,k) = sh_(:,k) + KG_s*inov(:,k);
    Ps = Ps_ - KG_s*Py*KG_s';
    s = Xh_s(:,k);
    
%%%  End of State Estimation  
    

    %Pw = Px; %eps*eye(size(Pw));
    %Rv_p = diag(max(0.9 * diag(Rv_p) , 1e-8));
    %Rv_s = diag(max(0.9 * diag(Rv_s) , 1e-8));
    %Rww_p = diag(max(0.9 * diag(Rww_p) , 1e-4));
    %Rww_s = diag(max(0.9 * diag(Rww_s) , 1e-10));
%%% End of Parameter Estimation
    
end

end

%  figure;plot(Xh_s(1,:))
% figure; plot(xh(2,:))
figure; plot(out_EEG,'k')
hold on,
plot((Xh_s(2,:) - Xh_s(3,:)),'r')

figure; plot(X(:,1),'k')
hold on,
plot(Xh_s(1,:),'r')

figure; plot(C2*ones(1,L),'k')
hold on,
plot(Xh_s(8,:),'r')



















