
function [theta_ex,alpha_ex,dthetaf,dalphaf,v_lqr,t] = PreData(AlphaFileName,ThetaFileName,sigma_theta,sigma_alpha)

% Loading data for training set
load(AlphaFileName,'alpha');
data_alpha=alpha(2,:);
alpha_ex=data_alpha'*pi/180;
t=alpha(1,:)';

load(ThetaFileName,'theta');
data_theta=theta(3,:);
theta_ex=data_theta'*pi/180;
theta_d=theta(2,:)'*pi/180;

dt=t(2);
%% V_lqr= k1(theta_d-theta)
k=[3.1623  -10.3958    0.7280   -0.3113];
v_lqr=k(1)*(theta_d-theta_ex);

%% Estimate Derivatives
% sigma_theta=0.05; % gain for theta signal
% sigma_alpha=0.05; % gain for alpha signal

est_theta = nleso(theta_ex',dt,sigma_theta);
thetaf = est_theta(1,:)';
dthetaf= est_theta(2,:)';

est_alpha = nleso(alpha_ex',dt,sigma_alpha);
alphaf = est_alpha(1,:)';
dalphaf= est_alpha(2,:)';

%%% plot
% figure()
% subplot(2,2,1)
% plot(t,theta_ex,t,thetaf)
% ylabel('theta')
% subplot(2,2,3)
% plot(t,dthetaf)
% ylabel('dtheta')
% 
% subplot(2,2,2)
% plot(t,alpha_ex,t,alphaf)
% ylabel('alpha')
% subplot(2,2,4)
% plot(t,dalphaf)
% ylabel('dalpha')
% print -depsc TestErrorRFJ.eps
end