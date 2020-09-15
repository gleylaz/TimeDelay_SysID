clear all;close all;clc

%% add libraries
addpath('Functions')

%% Generating Data
sigmaF=0.01;
sigmaX=0.01;

% time step setup
t_init=0;
T=5;
dt=0.0002;
t=transpose(t_init:dt:T);
nt=length(t);
tau=0.2;

% [initiLI position;initiLI velocity];
inits_tr=[1;0;0];
inits_tst=[-1;0.5;0];

% Training Data
Force_tr=10*sin(4*t)+40*sin(4*t.^2)+sigmaF*randn(nt,1);
[x_tr,v_tr]=sys(Force_tr,inits_tr,tau,t,dt);
x_tr=x_tr+sigmaX*randn(nt,1);%displacement

% Test Data
Force_tst=10*sin(2*t)+40*sin(8*t)+sigmaF*randn(nt,1);
[x_tst,v_tst]=sys(Force_tst,inits_tst,tau,t,dt);
x_tst=x_tst+sigmaX*randn(nt,1);%displacement

plot(t,x_tr,t,x_tst);
legend('Training Data','Test Data')
grid on
xlabel('$time$')
ylabel('$Displacement$')
title('Displacement Vs Time') 
print -depsc ResponseMassSpringDamper.eps

%% Algebraic derivatives with trapezoid integration

% Third Order
pf_tr=Al(3,2,Force_tr,t);
pf_tst=Al(3,2,Force_tst,t);

%% Shrinkage method joined with Bootsraping
MaxPol=5;
CrossedProducts=1;
CoulombFriction=1;

% make a vector of threshold values
numlambda = 50;
lambdastart = -6;
lambdaend = 0;
Lambda = logspace(lambdastart,lambdaend, numlambda);


K=10;
ratio=0.5;
L=round(ratio*nt);

parameters=cell(MaxPol,1);
avg_prms=cell(MaxPol,1);

for pol=1:MaxPol
    clear p_tr p_tst P_tr Pf_tr est MSE Estimated
    
    % Lib(x,dx,t,polyorder,CrossedProducts,CoulombFriction)
    p_tr=Lib(x_tr,v_tr,t,pol,CrossedProducts,CoulombFriction);
    p_tst=Lib(x_tst,v_tst,t,pol,CrossedProducts,CoulombFriction);
    
    for k=1:K
        idx=randperm(nt);
        P_tr=p_tr(idx(1:L),:);
        Pf_tr=pf_tr(idx(1:L),:);
        % Sparse Regression
        for i=1:numlambda
            est(:,i)=sparsifyDynamics(P_tr,Pf_tr,Lambda(i),1);
            MSE(i)=(p_tst*est(:,i)-pf_tst)'*(p_tst*est(:,i)-pf_tst)/nt;
        end

        % Sparsification lambda selection
        [MSE_min,Index]=min(MSE);
        lambdaMin(k) = Lambda(Index);
        Estimated(:,k)=sparsifyDynamics(P_tr,Pf_tr,lambdaMin(k),1);
    end
    
    avg_prms{pol}=mean(Estimated,2);
    AIC(pol)=(p_tst*avg_prms{pol}-pf_tst)'*(p_tst*avg_prms{pol}-pf_tst)/nt;
    parameters{pol}=Estimated;

end

%%
% figure()
% plot(Lambda,MSE,'LineWidth',1.5)
% xlabel('$\lambda$')
% xline(lambdaMin,'-.b','\lambda_{min}','LineWidth',1.5,'FontSize', 24);
% ylabel('$MSE$')
% grid on

%% Model Selection
semilogy(1:5,AIC)

%% feature selection
% logic_est=logical(Estimated);
% p_tau=p_tr(:,1);
% p_tr(:, logic_est== 0)= [];
% p=[p_tau,p_tr];
% 
% OLS=(p'*p)\(p'*Pf_tr)

%% Comparison coefficients

m=1;

k1=6;
k2=0;
k3=2;

c1=1;
c2=1;
c3=1;

k_p=2;
true_coef=[k_p*(tau^3)/6;
    m-k_p*0.5*(tau^2);
    k1-k_p;
    k3;
    c1+k_p*tau;
    c2;
    c3]

% comparison_table=table(true_coef,Estimated_OLS_Al)
