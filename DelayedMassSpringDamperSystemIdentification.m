
% Copyright 2020, All Rights Reserved
% Code by Ghazaale Leylaz
% For Paper, "A Robust Data-Driven Identification Algorithm for  
% Nonlinear Dynamical Systems with Time Delay"
% by Ghazaale Leylaz, Shuo Wang, and Jian-Qiao Sun

% A bootsraping sparse regeression algorithm with algebraic operation 
% to eleminate the effect of noise and estimatate signal derivatives 
% is developed for nonparametric indentification of a delayed nonlinear 
% mass-spring-damper system.


clear all;close all;clc

%% add libraries
addpath('Functions')
swEPSfigure
swFigSize

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

% % [initiLI position;initial velocity,initial acceleration];
% inits_tr=[1;0;0];
% inits_tst=[-1;0.5;0];
% 
% % Training Data
% Force_tr=10*sin(4*t)+40*sin(4*t.^2)+sigmaF*randn(nt,1);
% [x_tr,v_tr]=sys(Force_tr,inits_tr,tau,t,dt);
% x_tr=x_tr+sigmaX*randn(nt,1);%displacement
% 
% % Test Data
% Force_tst=10*sin(2*t)+40*sin(8*t)+sigmaF*randn(nt,1);
% [x_tst,v_tst]=sys(Force_tst,inits_tst,tau,t,dt);
% x_tst=x_tst+sigmaX*randn(nt,1);%displacement
% 
% plot(t,x_tr,t,x_tst);
% legend('Training Data','Test Data')
% grid on
% xlabel('$time$')
% ylabel('$Displacement$')
% title('Displacement Vs Time') 
% print -depsc ResponseMassSpringDamper.eps
% 
% %% Saving 
% Data.x_tr=x_tr;
% Data.v_tr=v_tr;
% Data.x_tst=x_tst;
% Data.v_tst=v_tst;
% Data.Force_tr=Force_tr;
% Data.Force_tst=Force_tst;
% Data.t=t;
% Data.nt=nt;
% save('DataMCK.mat','-struct','Data')

%% Loading system's responses and excitations data
load('DataMCK.mat')
%% Algebraic derivatives with trapezoid integration
% Third Order
pf_tr=Al(3,2,Force_tr,t);
pf_tst=Al(3,2,Force_tst,t);

%% Shrinkage method joined with Bootsraping
MaxPol=6;
CrossedProducts=1;
CoulombFriction=1;

% make a vector of threshold values
numlambda = 50;
lambdastart = -6;
lambdaend = 0;
Lambda = logspace(lambdastart,lambdaend, numlambda);

K=30; % 
ratio=0.5;
L=round(ratio*nt);

%% Algorithm
% parameters=cell(MaxPol,1);
% avg_prms=cell(MaxPol,1);
% std_prms=cell(MaxPol,1);
% 
% for pol=1:MaxPol
%     clear p_tr p_tst P_tr Pf_tr est MSE Estimated
%     
%     % Lib(x,dx,t,polyorder,CrossedProducts,CoulombFriction)
%     p_tr=Lib(x_tr,v_tr,t,pol,CrossedProducts,CoulombFriction);
%     p_tst=Lib(x_tst,v_tst,t,pol,CrossedProducts,CoulombFriction);
%     
%     for k=1:K
%         [P_tr,idx]=datasample(p_tr,L);
%         Pf_tr=pf_tr(idx,:);
%         % Sparse Regression
%         for i=1:numlambda
%             est(:,i)=sparsifyDynamics(P_tr,Pf_tr,Lambda(i),1);
%             MSE(i)=(p_tst*est(:,i)-pf_tst)'*(p_tst*est(:,i)-pf_tst)/nt;
%         end
% 
%         % Sparsification lambda selection
%         [MSE_min,Index]=min(MSE);
%         lambdaMin(k) = Lambda(Index);
%         Estimated(:,k)=sparsifyDynamics(P_tr,Pf_tr,lambdaMin(k),1);
%         
%         if Estimated(1,k)==0
%             logic_est=logical(Estimated(:,k));
%             p_tau=P_tr(:,1);
%             P_tr(:, logic_est== 0)= 0;
%             P_tr(:,1)=p_tau;
%             Estimated(:,k)=P_tr\Pf_tr;
%         end  
%     end
%     
%     avg_prms{pol}=mean(Estimated,2);
%     % Standard deviations
%     std_prms{pol}=std(Estimated,0,2);
%     Test_Error(pol)=(p_tst*avg_prms{pol}-pf_tst)'*(p_tst*avg_prms{pol}-pf_tst)/nt;
%     parameters{pol}=Estimated;
% 
% end
% 
% %% Save results
% Results.parameters=parameters;
% Results.avg_prms=avg_prms;
% Results.std_prms=std_prms;
% Results.Test_Error=Test_Error;
% save('ResultsMCK.mat','-struct','Results')

%% Load the saved results
load('ResultsMCK.mat')

%% Model Selection
[Test_Error_min,True_order]=min(Test_Error);
Estimated_prms=avg_prms{True_order}
Standard_deviation=std_prms{True_order};

%% Test Error plotting
figure
semilogy(1:MaxPol,Test_Error,'-*','LineWidth',1.5)
% xline(True_order,'-.b',['MSE=',num2str(Test_Error_min)],'LineWidth',1.5,'FontSize', 24);
ylabel('$MSE_{test}$')
xlabel('Polynomial Order')
hold on
plot(True_order,Test_Error_min,'or','LineWidth',3)
grid on
print -depsc TestErrorMCK.eps
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


