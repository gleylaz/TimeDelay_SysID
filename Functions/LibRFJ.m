
function [p] = LibRFJ(x,dx,t,polyorder,CrossedProducts,CoulombFriction)
 
%%% Algebraic Data Library

numTraj=size(x,2);

Pa=[];
Px=[];
Pv=[];
Pxv=[];
PCoulomb=[];

for n=1:numTraj
    % pa
    pa(:,1)=(Al(0,2,x(:,n),t)-4*Al(1,1,x(:,n),t)+2*Al(2,0,x(:,n),t));
    Pa=[Pa;pa];

    % px
    for i=1:polyorder
        px(:,i)=Al(2,2,x(:,n).^i,t);
    end
    Px=[Px;px];

    % pv
    pv(:,1)=Al(1,2,x(:,n),t)-2*Al(2,1,x(:,n),t);
    for j=2:polyorder
        pv(:,j)=Al(2,2,dx(:,n).^j,t);
    end
    Pv=[Pv;pv];

    % pxv
    for k=1:polyorder
        for l=1:polyorder
            pxv(:,(k-1)*polyorder+l)=Al(2,2,(x(:,n).^k).*(dx(:,n).^l),t);
        end
    end
    Pxv=[Pxv;pxv];

    % pCoulomb
    pCoulomb=Al(2,2,sign(dx(:,n)),t);
    PCoulomb=[PCoulomb;pCoulomb];
 
end



if CoulombFriction==0
    if CrossedProducts==1
        p=[Pa,Px,Pv,Pxv]; 
     elseif CrossedProducts==0
         p=[Pa,Px,Pv];
    end
    
elseif CoulombFriction==1
    if CrossedProducts==1
        p=[Pa,Px,Pv,Pxv,PCoulomb]; 
    elseif CrossedProducts==0
        p=[Pa,Px,Pv,PCoulomb];
    end
end

end


    