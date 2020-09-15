function [p] = Lib(x,dx,t,polyorder,CrossedProducts,CoulombFriction)
 
%%% Algebraic Data Library

p_tau=6*Al(2,0,x,t)-6*Al(1,1,x,t)+Al(0,2,x,t);

pa(:,1)=2*Al(3,0,x,t)-4*Al(2,1,x,t)+Al(1,2,x,t);

% px
for i=1:polyorder
    px(:,i)=Al(3,2,x.^i,t);
end

% pv
pv(:,1)=-2*Al(3,1,x,t)+Al(2,2,x,t);
for j=2:polyorder
    pv(:,j)=Al(3,2,dx.^j,t);
end

% pxv
for k=1:polyorder
    for l=1:polyorder
        pxv(:,(k-1)*polyorder+l)=Al(3,2,(x.^k).*(dx.^l),t);
    end
end

% pCoulomb
pCoulomb=Al(3,2,sign(dx),t);
    
if CoulombFriction==0
    if CrossedProducts==1
        p=[p_tau,pa,px,pv,pxv]; 
     elseif CrossedProducts==0
         p=[p_tau,pa,px,pv];
    end
    
elseif CoulombFriction==1
    if CrossedProducts==1
        p=[p_tau,pa,px,pv,pxv,pCoulomb]; 
    elseif CrossedProducts==0
        p=[p_tau,pa,px,pv,pCoulomb];
    end
end

end


    