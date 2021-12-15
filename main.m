
n=6;  % number of states
m=9; % number of measurements
p=3; % number of inputs
N=50; % number of instants;

phi=[1 0 0 0.01 0 0;
     0 1 0 0 0.01 0;
     0 0 1 0 0 0.01;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];
Gamma_u=[0.0000005 0 0;
          0 0.0000005 0;
          0 0 0.0000005;
          0.01 0 0;
          0 0.01 0;
          0 0 0.01]; 
Gamma_d=[0.01 0 0 0.0000005 0 0;
         0 0.01 0 0 0.0000005 0;
         0 0  0.01 0 0 0.0000005;
         0 0 0 0.01 0 0;
         0 0 0 0 0.01 0;
         0 0 0 0 0 0.01];
Cmat=[0 0 0.00003 0 0 0;
      0 0 0 0 0 0 ;
      0.00003 0 0 0 0 0;
      0 -9.81 0 0 0 0;
      9.81 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1];
Qmat=0.0001*eye(n);
Rmat=0.0001*eye(m);

%%%%%%%%%%%%%%%%%
%% True state generation
%%rand('state',0); % random number generator initialization
rng('default')
w=mvnrnd(zeros(N,n),Qmat); % process noise generation
w=w'; % column corresponds to time instants
v=mvnrnd(zeros(N,m),Rmat); % measurement noise generation
v=v'; 
xt=zeros(n,N); 
xt(:,1)=0.05*ones(n,1); % the initial true state
y=zeros(m,N);
u=zeros(p,N);

for k=1:N
   u(1,k)=0.5*exp(-0.01*k)-0.5*exp(-0.005*k);
   u(2,k)=-0.5*exp(-0.01*k)-0.5*exp(-0.005*k);
   u(3,k)=-0.5*exp(-0.005*k);
end

for k=1:N-1
    xt(:,k+1)=phi*xt(:,k) + Gamma_u*u(:,k)+ w(:,k);
    y(:,k+1)=Cmat*xt(:,k+1) + v(:,k+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% The Kalman Filter
kxhat_est=zeros(n,N);
kxhat_est(:,1)= -0.05*ones(n,1); % initializing xhat(0|0)
kPest=0.0001*eye(n);% initializing P(0|0)
kP=zeros(n,n);
kP(:,:,1)=kPest;
kPpred=0.0001*eye(n); % initializing P(0|-1) just for completeness. It is not required in computations 
ktracePest(1)=trace(kPest);
kspecPest(1)=max(eig(kPest));
ktracePpred(1)=trace(kPpred);
kspecPpred(1)=max(eig(kPpred));

for k=1:N-1
    kxhat_pred(:,k+1)=phi*kxhat_est(:,k)+Gamma_u*u(:,k);
    kPpred=phi*kPest*phi' + Qmat;
    ktracePpred(k+1)=trace(kPpred);
    kspecPpred(k+1)=max(eig(kPpred));
    kL=kPpred*Cmat'*inv(Cmat*kPpred*Cmat'+Rmat);
    kPest=(eye(n)-kL*Cmat)*kPpred;
    kP(:,:,k+1)=kPest;
    ke(:,k+1)=y(:,k+1)-Cmat*kxhat_pred(:,k+1);
    kxhat_est(:,k+1)=kxhat_pred(:,k+1)+kL*ke(:,k+1);
    ktracePest(k+1)=trace(kPest);
    kspecPest(k+1)=max(eig(kPest));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extended Kalman Filter
exhat_est=zeros(n,N);
exhat_est(:,1)= -0.05*ones(n,1); % initializing xhat(0|0)
ePest=0.0001*eye(n); % initializing P(0|0)
eP=zeros(n,n);
eP(:,:,1)=ePest;
ePpred=0.0000001*eye(n); % initializing P(0|-1) just for completeness. It is not required in computations 
etracePest(1)=trace(ePest);
especPest(1)=max(eig(ePest));
etracePpred(1)=trace(ePpred);
especPpred(1)=max(eig(ePpred));

for k=1:N-1
    % Use kxhat_pred(:,k+1)=kA(exhat_est(:,k))*exhat_est(:,k)+Gamma_u*u(:,k);
    exhat_pred(:,k+1)=Areal(exhat_est(:,k))*exhat_est(:,k)+Gamma_u*u(:,k);
    % c2d is used find phi and gamma_d for this filter
     A = kA(exhat_est(:,k));
     B = [0 0 0;
          0 0 0;
          0 0 0;
          1 0 0;
          0 1 0;
          0 0 1];
     F=eye(n);
     C = [0 0 0 0 0 0];
     D = 0;
     % The disturbance matrix is combined with input
     % for using c2d
    sys1 = ss(A,[B F],C,D);
    sysd1 = c2d(sys1,0.1,'tustin');
    [ephi,big,c,d]=ssdata(sysd1);
    kGamma_d = big(:,4:end) ;
    ePpred=ephi*ePest*ephi' + kGamma_d*Qmat*kGamma_d';
    etracePpred(k+1)=trace(ePpred);
    especPpred(k+1)=max(eig(ePpred));
    % Compute Cmat using jacobian of G()
    kCmat=kCm(exhat_pred(:,k+1));
    eL=ePpred*kCmat'*inv(kCmat*ePpred*kCmat'+ Rmat );
    ePest=(eye(n)-eL*kCmat)*ePpred;
    eP(:,:,k+1)=ePest;
    ee(:,k+1)=y(:,k+1)-kG(exhat_pred(:,k+1));
    exhat_est(:,k+1)=exhat_pred(:,k+1)+eL*ee(:,k+1);
    etracePest(k+1)=trace(ePest);
    especPest(k+1)=max(eig(ePest));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Unscented Kalman Filter
M=21; % 21=6(number of states)+6(size of disturbance)+9(size of output)
uxhat_est=zeros(n,N);
uxhat_est(:,1)= -0.05*ones(n,1);% initializing xhat(0|0)
uxhat_pred=zeros(n,N);
uPest=0.0001*eye(n) ;% initializing P(0|0)
uP=zeros(n,n);
uP(:,:,1)=ePest;
uPpred=0.0001*eye(M); % initializing P(0|-1) just for completeness. It is not required in computations ()
utracePest(1)=trace(uPest);
uspecPest(1)=max(eig(uPest));
utracePpred(1)=trace(uPpred);
uspecPpred(1)=max(eig(uPpred));
% Defining parameters
kappa=1;
r=21+kappa;
rho=sqrt(r);
weights=(1/(2*r))*ones(1,2*M+1);
weights(1)=kappa/r;
% Estimation samples 
p=eye(M);
uxhat_est_samples=zeros(21,43);
uxhat2_est_samples=zeros(21,43);
uy_samples=zeros(m,13);

Pespe=zeros(n,m);
Pee=zeros(m,m);
uxhat_pred_samples=zeros(6,43);

for k=1:N-1
    uxhat_est_samples(:,1)=[uxhat_est(:,k);zeros(n,1);zeros(m,1)];
    % Augmented P matrix
    uPest=[uPest,zeros(6,6),zeros(6,9);
        zeros(6,6), Qmat, zeros(6,9);
        zeros(9,6), zeros(9,6), Rmat];
    for i=1:21
    uxhat_est_samples(:,i+1)=uxhat_est_samples(:,1)+rho*sqrtm(uPest)*p(:,i);
    uxhat_est_samples(:,i+22)=uxhat_est_samples(:,1)-rho*sqrtm(uPest)*p(:,i);
    end
    % Sample propogation
    for i=1:43
        uxhat_pred_samples(:,i)=Areal(uxhat_est_samples(1:n,i))*uxhat_est_samples(1:n,i)+Gamma_u*u(:,k)+Gamma_d*uxhat_est_samples(n+1:2*n,i);
    end 
    % xhat_pred
    for i=1:43
        uxhat_pred(:,k+1)=uxhat_pred(:,k+1)+weights(i)*uxhat_pred_samples(:,i);
    end
    Pepseps=zeros(n,n);
    for i=1:43
       Pepseps=Pepseps+weights(i)*(uxhat_pred_samples(1:n,i)-uxhat_pred(:,k+1))*transpose(uxhat_pred_samples(1:n,i)-uxhat_pred(:,k+1));
    end
    % Pepseps treated as P(k+|k)
    utracePpred(k+1)=trace(Pepseps);
    uspecPpred(k+1)=max(eig(Pepseps));
    % Output Samples
    for i=1:43
        uy_samples(:,i)=kG(uxhat_pred_samples(1:n,i))+v(:,k+1);
    end
    % Output Sample Mean
    uy_samples_mean=zeros(9,1);
    for i=1:43
        uy_samples_mean=uy_samples_mean+weights(i)*uy_samples(:,i);
    end
    % P-epsilon-e
    for i =1:43
        Pespe=Pespe+weights(i)*(uxhat_pred_samples(1:n,i)-uxhat_pred_samples(1:n,1))*(uy_samples(:,i)-uy_samples_mean)';
    end
    % P-e-e
    for i =1:43
        Pee=Pee+weights(i)*(uy_samples(:,i)-uy_samples_mean)*(uy_samples(:,i)-uy_samples_mean)';
    end
    % Update step
    uL=Pespe/Pee;
    uPest=Pepseps-uL*Pee*transpose(uL);
    uP(:,:,k+1)=uPest;
    ue(:,k+1)=y(:,k+1)-uy_samples_mean;
    uxhat_est(:,k+1)=uxhat_pred(:,k+1)+uL*ue(:,k+1);
    utracePest(k+1)=trace(uPest);
    uspecPest(k+1)=max(eig(uPest));
end
 esterror1=xt(:,1:20)-kxhat_est(:,1:20);
 esterror2=xt(:,1:20)-exhat_est(:,1:20);
 esterror3=xt(:,1:20)-uxhat_est(:,1:20);
% b1=[3*(0.01),bound(esterror1(6,1:19))];
% b2=[3*(0.01),bound(esterror2(6,1:19))];
% b3=[3*(0.01),bound(esterror3(6,1:19))];
plot([1:N],exhat_est(3,:),[1:N],xt(3,:) );
title('Sensitivity Analysis')
xlabel('Iterations')
ylabel('state')
legend('UKF','True state')
emean1=mean(ke,2);
emean2=mean(ee,2);
emean3=mean(ue,2);

covk1=cov(ke(1,:));
covk2=cov(ke(2,:));
covk3=cov(ke(3,:));
covk4=cov(ke(4,:));
covk5=cov(ke(5,:));
covk6=cov(ke(6,:));
covk7=cov(ke(7,:));
covk8=cov(ke(8,:));
covk9=cov(ke(9,:));

cove1=cov(ee(1,:));
cove2=cov(ee(2,:));
cove3=cov(ee(3,:));
cove4=cov(ee(4,:));
cove5=cov(ee(5,:));
cove6=cov(ee(6,:));
cove7=cov(ee(7,:));
cove8=cov(ee(8,:));
cove9=cov(ee(9,:));


covu1=cov(ue(1,:));
covu2=cov(ue(2,:));
covu3=cov(ue(3,:));
covu4=cov(ue(4,:));
covu5=cov(ue(5,:));
covu6=cov(ue(6,:));
covu7=cov(ue(7,:));
covu8=cov(ue(8,:));
covu9=cov(ue(9,:));


kbeta=zeros(1,N);
ebeta=zeros(1,N);
ubeta=zeros(1,N);

for i=1:N
    kbeta(i)=(kxhat_est(:,i)-mean(kxhat_est,2))'*inv(kP(:,:,i))*(kxhat_est(:,i)-mean(kxhat_est,2));
end

for i=1:N
    ebeta(i)=(exhat_est(:,i)-mean(exhat_est,2))'*inv(eP(:,:,i))*(exhat_est(:,i)-mean(exhat_est,2));
end

for i=1:N
    ubeta(i)=(uxhat_est(:,i)-mean(uxhat_est,2))'*inv(kP(:,:,i))*(uxhat_est(:,i)-mean(uxhat_est,2));
end
alpha=0.05;
kai1 = chi2inv(alpha,n);
kai2 = chi2inv(1-alpha,n);



%   plot([1:N],ubeta,[1:N], kai1*ones(1,N),[1:N], kai2*ones(1,N));
%   title('Normalised Estimation Error Squared for ukf')
%   xlabel('Iterations')
%   ylabel('Beta')
%legend('Estimation Error-kf','Upper limit - kf','Lower limit - kf','Estimation Error-ekf', 'Upper limit - ekf','Lower limit - ekf', 'Estimation Error-ukf', 'Upper limit - ukf','Lower limit - ukf')
