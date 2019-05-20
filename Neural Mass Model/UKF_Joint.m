

T = 0.01;
N=length(0:T:200);
u = zeros(1,N);
u(20/T:100/T) = 0.1;
B = [1 0 0 0]';
C = B';


Xh_s=[];
Xh_=[];

L2=1;%size(observation);
Obs_dim=L2;
s_dim = 7;
nsp_s = 2*s_dim + 1;

%%% Weighths
alpha = 1;
ktau = 0;%3-L1;
beta = 3;
lambda_s = alpha^2 *(s_dim+ktau) - s_dim;
W_s = [lambda_s 0.5 0]/(s_dim+lambda_s); 
W_s(3) = W_s(1) + (1-alpha^2) + beta;
gamma_s = sqrt(s_dim+lambda_s);
%%%%%%%%%%
Rv_s = 1e-10*eye(s_dim);
Rww_s = 10*eye(Obs_dim);
Sv_s = chol(Rv_s)';
%Sw = chol(Rww_s)';
v=-64.9964; m=0.0530; h=0.5960; n=0.3177;
s = [v;n;m;h;gNA_main+2;gK_main-2;gL_main];
Xh_s(:,1) = s;
Ps = 0.01*eye(size(s,1));

aa=0.9;
for j=1:1
for k=2:N
 
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
        aaa = F_value_dual(Z(:,i),cM_main, Z(5:7,i), ENA, EK, EL, T) + T*B*(u(k)/A_M * 1e-6)/(cM_main*1e-2); % F(X_sigma(:,i),input);
        X_(:,i)=[aaa;Z(5:7,i)];% + [sqrt(1)*randn(1);0;0;0];
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
    Y_ = X_(1,:); % observation function
    yh_(:,k) = W_s(1)*Y_(:,1) + W_s(2)*sum(Y_(:,2:nsp_s),2);
    temp2 = Y_(1,:) - cvecrep(yh_(1,k),nsp_s);
    Py  = W_s(3)*temp2(:,1)*temp2(:,1)' + W_s(2)*temp2(:,2:nsp_s)*temp2(:,2:nsp_s)' + Rww_s;
   
    Psy = W_s(3)*temp1(:,1)*temp2(:,1)' + W_s(2)*temp1(:,2:nsp_s)*temp2(:,2:nsp_s)';
    KG_s = Psy / Py;
    inov(:,k) = y(:,k) - yh_(1,k);
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

 figure;plot(Xh_s(1,:))
% figure; plot(xh(2,:))























