

Xh_s=[];
Xh_=[];

obs = out_EEG;
L2=1;%size(observation);
Obs_dim=L2;
s_dim = 8; % x1:x6 & C1:C6
nsp_s = 2*s_dim;% + 1;

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

%Sw = chol(Rww_s)';
x1=0; x2=0; x3=0; x4=0; x5=0; x6=0; c1 = 400; c2 = 20;
s = [x1;x2;x3;x4;x5;x6;c1;c2];
k_init = 1000;
Xh_s(:,k_init) = s;
%% Tuning the Noise variances
st1 = std(X(5100:5200,1)); st2 = std(X(5100:5200,2)); st3 = std(X(5100:5200,3));
st4 = std(X(5100:5200,4)); st5 = std(X(5100:5200,5)); st6 = std(X(5100:5200,6));
Rv_s = diag([1e-4*[st1,st2,st3,st4,st5,st6],1e-2*c1,1e-2*c2]);%diag([1e-8*ones(1,6),1e-2*c1,1e-2*c2]);%5e-10*eye(s_dim);
Rww_s = 1e-2*eye(Obs_dim);
Sv_s = chol(Rv_s)';
Ps = 1e-2*diag([1e-4*[st1,st2,st3,st4,st5,st6],1e-2*c1,1e-2*c2]);%diag([1e-6*ones(1,6),1e-6*c1,1e-6*c2]);%1*eye(size(s,1));
%%
Ds = 6; % Number of differential equations describing model, also the number of fast states to be estiamted

Dp = 2; % Number of parameters to be estimated, also refered to as slow states

Dk = 0; %If set to 1 the mean of the stochastic input will be estimated % Note that if Input_mean_variation is not zero than Dk should be set to one to allow tracking of the input mean

Dy =1; % Number of observable outputs from the simulation

Dx = Ds+Dp+Dk; % Number of dimensions of augmented state matrix, Note that estimated parameters and inputs are now considered to be 'slow states' in the estimation procedure

kappa =0; % Varibale used to define the relative contribution of the mean on the propogation of states, and adjustment of the variance of the sigma points drawn from the Gaussian distribution
if kappa > 0
    W_s = [kappa/(Dx+kappa) ones(1,2*s_dim)*1/(2*(kappa+Dx))];
    Number_sigma = 2*Dx+1;
else
    W_s = ones(1,2*s_dim)*1/(2*(kappa+Dx));
    Number_sigma = 2*Dx;
end
%%
aa=0.9;
for j=1:1
for k=k_init+1:length(obs)
 
%Sw = chol(Rww)';   
%%%    State Estimation
%     Z=[];
%     Sx_s = chol(Ps)';
%     Z    = cvecrep(s, nsp_s);                                                   
%     sSz  = gamma_s * Sx_s;
%     sSzs = [sSz -sSz];
%     Z(:,2:nsp_s) = Z(:,2:nsp_s) + sSzs;           % build sigma-point set
[Sigma(:,:,k) err] = Unscented_transform(Dx,Ps(:,:),Xh_s(:,k-1),kappa);
Z = Sigma(:,:,k);
    X_=[];
    Y_=[];
    temp1=[];
    for i=1:nsp_s
        p.C1 = Z(7,i); p.C2 = Z(8,i); p.C3 = C3; p.C4 = C4;
        aaa = JR_Model(Z(1:Ds,i),mu(k-1),Ip(k-1),Ii(k-1),p); % F(X_sigma(:,i),input);
        X_(:,i)=[aaa;Z(Ds+1:Ds+Dp,i)] + 0*randn(s_dim,1);
    end
    sh_(:,k) = sum(bsxfun(@times,X_',W_s'))';%W_s(1)*X_(:,1) + W_s(2)*sum(X_(:,2:nsp_s),2);% + 0.01*randn(size(xh_(:,k-1)));
    temp1 = X_ - cvecrep(sh_(:,k),nsp_s);
    Ps_ = bsxfun(@times,temp1',W_s')'*temp1' + Rv_s;%W_s(3)*temp1(:,1)*temp1(:,1)' + W_s(2)*temp1(:,2:nsp_s)*temp1(:,2:nsp_s)'+ Rv_s;
    %%% augmented sigma points
%     Z=[];
%     Sx_s = chol(Ps_)';
%     Z    = cvecrep(sh_(:,k), nsp_s);                                                   
%     sSz  = gamma_s * Sx_s;
%     sSzs = [sSz -sSz];
%     Z(:,2:nsp_s) = Z(:,2:nsp_s) + sSzs;
    %%%
    Y_=[];
    %yh_=[];
    temp2=[];
    Y_ = X_(2,:) - X_(3,:); %Y_ = Y_ - mean(Y_);%X_(2,:) - X_(3,:);%  % observation function
    yh_(:,k) = sum(W_s.*Y_);%W_s(1)*Y_(:,1) + W_s(2)*sum(Y_(:,2:nsp_s),2);
    temp2 = Y_(1,:) - cvecrep(yh_(1,k),nsp_s);
    Py = bsxfun(@times,temp2',W_s')'*temp2'+ Rww_s; %W_s(3)*temp2(:,1)*temp2(:,1)' + W_s(2)*temp2(:,2:nsp_s)*temp2(:,2:nsp_s)' + Rww_s;
   
    Psy = bsxfun(@times,temp1',W_s')'*temp2';%W_s(3)*temp1(:,1)*temp2(:,1)' + W_s(2)*temp1(:,2:nsp_s)*temp2(:,2:nsp_s)';
    KG_s = Psy* Py ^-1;
    inov(:,k) = obs(k) - yh_(:,k);
    Xh_s(:,k) = sh_(:,k) + KG_s*inov(:,k);
    Ps = Ps_ - KG_s*Py*KG_s';
    s = Xh_s(:,k);
    
    %Pw = Px; %eps*eye(size(Pw));
    %Rv_p = diag(max(0.9 * diag(Rv_p) , 1e-8));
    %Rv_s = diag(max(0.9 * diag(Rv_s) , 1e-8));
    %Rww_p = diag(max(0.9 * diag(Rww_p) , 1e-4));
    %Rww_s = diag(max(0.9 * diag(Rww_s) , 1e-10));
%%% End of Parameter Estimation
    
end

end
%% Figures
%  figure;plot(Xh_s(1,:))
% figure; plot(xh(2,:))
indx = 4000:16000;
tt2 = indx(1)*dt:dt:(indx(end))*dt;

figure; plot(tt2,out_EEG(indx),'k')
hold on,
out_Est = (Xh_s(2,indx) - Xh_s(3,indx));
plot(tt2,out_Est,'r') % plot( yh_(indx),'r')
xlabel('Time (Sec)')
ylabel('Voltage (mV)')
title('EEG')
legend('Recorded','Estimated')

figure; plot(tt2,X(indx,1),'k')
hold on,
plot(tt2,Xh_s(1,indx),'r')
xlabel('Time (Sec)')
ylabel('Amp')
title('Hidden State(1)')
legend('Original','Estimated')

figure; plot(tt,C1*ones(1,length(tt)),'k')
hold on,
plot(tt,Xh_s(7,:),'r')
title('Parameter C1')
xlabel('Time (Sec)')
ylabel('Amp')
legend('Original','Estimated')

figure; plot(tt,C2*ones(1,length(tt)),'k')
hold on,
plot(tt,Xh_s(8,:),'r')
title('Parameter C2')
xlabel('Time (Sec)')
ylabel('Amp')
legend('Original','Estimated')













