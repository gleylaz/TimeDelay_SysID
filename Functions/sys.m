
% The function generates data 
function [x,v]=sys(Force,inits,tau,t,dt)

% System paratameters for modeLIng the system: mddx=-k1*x-k2*x^2-k3*x^3-c1*dx-c2*x^2*dx-c3*dx^3
m=1;

k1=6;
k2=0;
k3=2;

c1=1;
c2=1;
c3=1;

k_p=2;


nt=length(t);
n_delay=round(tau/dt);

% inits=[initiLI position;initiLI velocity;initiLI velocity];
x(:,1)=inits(1:2)';
dx(:,1)=inits(2:3)';

% Euler method to solve the motion equation
for n=2:n_delay
    dx(:,n)=[x(2,n-1);(Force(n-1)-k1*x(1,n-1)^1-k2*(x(1,n-1)^2)-k3*(x(1,n-1)^3)-c1*x(2,n-1)-c2*x(2,n-1)*(x(1,n-1)^2)-c3*(x(2,n-1)^3))/m];
    x(:,n)=x(:,n-1)+dx(:,n)*dt;
end

for n=n_delay+1:nt
    dx(:,n)=[x(2,n-1);(Force(n-1)-k1*x(1,n-1)^1-k2*(x(1,n-1)^2)-k3*(x(1,n-1)^3)-c1*x(2,n-1)-c2*x(2,n-1)*(x(1,n-1)^2)-c3*(x(2,n-1)^3)+k_p*x(1,n-n_delay))/m];
    x(:,n)=x(:,n-1)+dx(:,n)*dt;
end

v=x(2,:)';% speed
x=x(1,:)';% displacement

end