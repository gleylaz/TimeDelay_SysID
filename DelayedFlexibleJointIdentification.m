
% Copyright 2020, All Rights Reserved
% Code by Ghazaale Leylaz
% For Paper, "A Robust Data-Driven Identification Algorithm for  
% Nonlinear Dynamical Systems with Time Delay"
% by Ghazaale Leylaz, Shuo Wang, and Jian-Qiao Sun

% A bootsraping sparse regeression algorithm with algebraic operation 
% to eleminate the effect of noise and estimatate signal derivatives 
% is developed for nonparametric indentification of a delayed rotary 
% flexible joint under a pre-designed LQR controller.

clear all;close all;clc

%% import functions and data folders
addpath('Functions')
addpath('Data')

swEPSfigure
swFigSize

%% Loading data for train set
[theta_ex,alpha_ex,dthetaf,dalphaf,v_lqr,t] = ...
    PreData('alpha_Sin_2_0.66.mat','theta_Sin_2_0.66.mat',0.05,0.05);

%% Loading data for test set
[theta_ex_test,alpha_ex_test,dthetaf_test,dalphaf_test,v_lqr_test,t_test] = ...
    PreData('alpha_Square_1_0.66.mat','theta_Square_1_0.66.mat',0.05,0.05);

%% Plotting responses
%%% plot train signals
figure()
subplot(3,2,[1,2])
plot(t,v_lqr,'r')
ylabel('$V_{m}$')
subplot(3,2,3)
plot(t,theta_ex)
ylabel('$\theta$')
subplot(3,2,5)
plot(t,dthetaf)
xlabel('$t (s)$')
ylabel('$\dot{\theta}$')

subplot(3,2,4)
plot(t,alpha_ex)
ylabel('$\alpha$')
subplot(3,2,6)
plot(t,dalphaf)
ylabel('$\dot{\alpha}$')
xlabel('$t (s)$')
print -depsc TrainSignalsRFJ.eps

%%% plot test signals
figure()
subplot(3,2,[1,2])
plot(t,v_lqr_test,'r')
ylabel('$V_{m}$')
subplot(3,2,3)
plot(t,theta_ex_test)
ylabel('$\theta$')
subplot(3,2,5)
plot(t,dthetaf_test)
xlabel('$t (s)$')
ylabel('$\dot{\theta}$')

subplot(3,2,4)
plot(t,alpha_ex_test)
ylabel('$\alpha$')
subplot(3,2,6)
plot(t,dalphaf_test)
ylabel('$\dot{\alpha}$')
xlabel('$t (s)$')
print -depsc TestSignalsRFJ.eps

%% True coefficients and system parameters

% time delay
dt=t(2)-t(1);
nt=length(t);
tau=200*dt;

m_1=0.064;%kg
m_2=0.030;%kg
L_1=0.298;%m
L_2=0.156;%m
d=0.235;%m
J_l=((m_1*(L_1^2))/3)+((m_2*(L_2^2))/12)+(m_2*(d^2));

K_s=1.3; % gained from natural frequency of the system (experimental value)

% Servo constants (from its manual)
B_eq=0.015;%high-gear equivalent viscous damping coefficient
%B_eq=0.025;

B_c=0;
B_l=0;

%J_eq=0.0021;%with load
J_eq=0.002084179192918;% precise value from the devise itself
J_s=K_s*(J_l+J_eq)/J_l;

n_g=0.9; % Gear Efficiency

k_g=70; % High-gear total gear ratio
n_m=0.69; % Motor Efficiency

k_t=7.68e-3;
k_m=7.68e-3;

R_m=2.6; %ohm

a=n_g*k_g*n_m*k_t/R_m;
b=-n_g*(k_g^2)*n_m*k_t*k_m/R_m;

K_LQR=[3.1623  -10.3958    0.7280   -0.3113]; % LQR controller gains

%% Algebraic identification method

% Excitation
bf=a.*Al(2,2,v_lqr,t)+b.*(Al(1,2,theta_ex,t)-2*Al(2,1,theta_ex,t));
bf_test=a.*Al(2,2,v_lqr_test,t)+b.*(Al(1,2,theta_ex_test,t)-2*Al(2,1,theta_ex_test,t));

% Theta Equation
CrossedProducts=0;
CoulombFriction_theta=1;
CoulombFriction_alpha=0;

%% Settings
% Lambda grid for sparse regression
numlambda = 100;
lambdastart = -10;
lambdaend = 0;
Lambda = logspace(lambdastart,lambdaend, numlambda);

% Sparse Regression with Bootsraping
MaxPol=7;

K=10; % generate K data sample from the orginal data set 
ratio=0.5;
L=round(ratio*nt);

%%
% parameters1=cell(MaxPol,1);
% parameters2=cell(MaxPol,1);
% 
% avg_prms1=cell(MaxPol,1);
% avg_prms2=cell(MaxPol,1);
% 
% std_prms1=cell(MaxPol,1);
% std_prms2=cell(MaxPol,1);
% 
% for pol=5:MaxPol
%     
%     clear p_theta p_alpha p_theta_test p_alpha_test p1 p1_test ...
%         p2 p2_test P1 P2 est1 est2 MSE Estimated1 Estimated2 
%     
%     p_theta=LibRFJ(theta_ex,dthetaf,t,pol,CrossedProducts,CoulombFriction_theta);
%     p_alpha=LibRFJ(alpha_ex,dalphaf,t,pol,CrossedProducts,CoulombFriction_alpha);
% 
%     p_theta_test=LibRFJ(theta_ex_test,dthetaf_test,t,pol,CrossedProducts,CoulombFriction_theta);
%     p_alpha_test=LibRFJ(alpha_ex_test,dalphaf_test,t,pol,CrossedProducts,CoulombFriction_alpha);
% 
%     % Theta Equation 
%     p1=[p_theta(:,1:end-1),p_alpha(:,2:end),p_theta(:,end)];
%     p1_test=[p_theta_test(:,1:end-1),p_alpha_test(:,2:end),p_theta_test(:,end)];
% 
%     % Alpha Equation
%     p2=[p_alpha(:,1),p_theta(:,1:end-1),p_alpha(:,2:end),p_theta(:,end)];
%     p2_test=[p_alpha_test(:,1),p_theta_test(:,1:end-1),p_alpha_test(:,2:end),p_theta_test(:,end)];
% 
%     for k=1:K
%         [P1,idx]=datasample(p1,L);
%         P2=p2(idx,:);
%         Pf=bf(idx,:);
%         
%         % Sparse Regression
%         for i=1:numlambda
%             est1(:,i)=sparsifyDynamics(P1,Pf,Lambda(i),1);
%             MSE1(i)=(p1_test*est1(:,i)-bf_test)'*(p1_test*est1(:,i)-bf_test)/nt;
%             est2(:,i)=sparsifyDynamics(P2,-Pf,Lambda(i),1);
%             MSE2(i)=(p2_test*est2(:,i)+bf_test)'*(p2_test*est2(:,i)+bf_test)/nt;
%         end
% 
%         % Model Selection
%         [MSE_min1,Index1]=min(MSE1);
%         lambdaMin1(k) = Lambda(Index1);
%         Estimated1(:,k)=sparsifyDynamics(P1,Pf,lambdaMin1(k),1);
%         
%         [MSE_min2,Index2]=min(MSE2);
%         lambdaMin2(k) = Lambda(Index2);
%         Estimated2(:,k)=sparsifyDynamics(P2,-Pf,lambdaMin1(k),1);
%         
%         if Estimated2(2,k)==0 % to make a constrain to keep the time delay term 
%             logic_est=logical(Estimated2(:,k));
%             p_tau=P2(:,2);
%             P2(:, logic_est== 0)= 0;
%             P2(:,2)=p_tau;
%             Estimated2(:,k)=P2\-Pf;
%         end  
%     end
%    
%     % mean values
%     avg_prms1{pol}=mean(Estimated1,2);
%     avg_prms2{pol}=mean(Estimated2,2);
%     
%     % Standard deviations
%     std_prms1{pol}=std(Estimated1,0,2);
%     std_prms2{pol}=std(Estimated2,0,2);
%     
%     Test_Error1(pol)=(p1_test*avg_prms1{pol}-bf_test)'*(p1_test*avg_prms1{pol}-bf_test)/nt;
%     Test_Error2(pol)=(p2_test*avg_prms2{pol}+bf_test)'*(p2_test*avg_prms2{pol}+bf_test)/nt;
%     
%     parameters1{pol}=Estimated1;
%     parameters2{pol}=Estimated2;
% end
% 
% %% Save results
% Results.parameters1=parameters1;
% Results.parameters2=parameters2;
% 
% Results.avg_prms1=avg_prms1;
% Results.avg_prms2=avg_prms2;
% 
% Results.std_prms1=std_prms1;
% Results.std_prms2=std_prms2;
% 
% Results.Test_Error1=Test_Error1;
% Results.Test_Error2=Test_Error2;
% 
% save('ResultsRFL.mat','-struct','Results')

%% Load the saved results
load('ResultsRFL.mat')

%% Model selection
[Test_Error_min1,True_order_1]=min(Test_Error1);
[Test_Error_min2,True_order_2]=min(Test_Error2);

Estimated_theta=avg_prms1{True_order_1};
Estimated_alpha=avg_prms2{True_order_2};

%% Test Error plotting
figure

subplot(2,1,1)
semilogy(1:MaxPol,Test_Error1,'*','LineWidth',1.5)
% xline(True_order_1,'-.b','LineWidth',1.5,'FontSize', 24);
ylabel('$MSE_{test,\alpha}$')
hold on
plot(True_order_1,Test_Error_min1,'or','LineWidth',3)
grid on

subplot(2,1,2)
semilogy(1:MaxPol,Test_Error2,'*','LineWidth',1.5)
% xline(True_order_2,'-.b','LineWidth',1.5,'FontSize', 24);
ylabel('$MSE_{test,\theta}$')
xlabel('Polynomial Order')
hold on
plot(True_order_2,Test_Error_min2,'or','LineWidth',3)
grid on

print -depsc TestErrorRFJ.eps
%%
Excat_Value1=[J_eq+0.5*a*K_LQR(1)*(tau^2);
              B_eq-a*K_LQR(1)*tau;
              -K_s;
              B_c];
          
Excat_Value2=[J_eq;
              -0.5*a*K_LQR(1)*(tau^2);
              a*K_LQR(1)*tau-B_eq;
              K_s*(J_l+J_eq)/J_l;
              B_c];
%% 
Estimated_tau=sqrt(-2*(Estimated_alpha(2))/(a*K_LQR(1)))
