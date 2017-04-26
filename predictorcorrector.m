clc;
clear all;

%% Input Data.....

rr=1;              % alpha.....
tau=-0.2;           % scaling for continuation.......
scale=1;           % scaling for Ql.........

%% Input data got from first freeing only Pl as free parameter.....

V2=0.5;                % V2........
alpha=6.8661e-04;      % Alpha........
Lambda=-1.0034;        % Lambda 1.....

%% 

hat_Lambda=1;      % Lambda 2.....


Tolee = 1;  
Iterr = 1;
counterr = 0;
while (Tolee > 1e-4)

%% Input Data for the first point on the boundary..........


z1=[alpha;V2;Lambda;hat_Lambda];

%% Power balance........

Pl=0.1;
Ql=(-0.25);

Ldfactor=1;


%% Predictor step taking z1 values ..........

% Argumenting extended Jacobian with an additional row and column for Lambda 2

Ext_JacobP=[V2*cos(alpha) sin(alpha) Pl*Ldfactor 0; V2*sin(alpha)  2*V2-cos(alpha) 0 -Ql; -2*(V2^2)*sin(alpha) (4*V2*cos(alpha))-1 0 0; 0 0 0 1];

right_vector=[0 0 0 1];       % Right hand side

tang_vector=Ext_JacobP/right_vector;        % Tangent Vector........

alpha=alpha+tau*tang_vector(1,1);
V2=V2+tau*tang_vector(2,1);
Lambda=Lambda+tau*tang_vector(3,1);
hat_Lambda=hat_Lambda+tau*tang_vector(4,1);


%zp=[0.0007;0.4875;-1.0034;0.9500];


zxx=[alpha;V2;Lambda;hat_Lambda];


%zp=z1+tau*tang_vector;                      % Predicted step on the boundry ..........


%% Corrector...................

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-6)
    
    ff1=((V2*sin(alpha)))+Pl*(1+Ldfactor*Lambda);
    ff2=(((V2^2) -V2*cos(alpha)))- hat_Lambda*Ql*(1);
    ff3=2*V2^2 *cos(alpha) -V2;
    ff4=(zxx-z1)'*tang_vector-tau;

    
    %norm(zxx-z1)^2-0.06;
        
    Ext_Jacob=[V2*cos(alpha) sin(alpha) Pl*Ldfactor 0; V2*sin(alpha)  2*V2-cos(alpha) 0 -Ql; -2*(V2^2)*sin(alpha) (4*V2*cos(alpha))-1 0 0;2*(zxx(1,1)-tang_vector(1,1)) 2*(zxx(2,1)-tang_vector(2,1)) 2*(zxx(3,1)-tang_vector(3,1)) 2*(zxx(4,1)-tang_vector(4,1))];
   
    
    dM=[ff1;ff2;ff3;ff4];
    
    corr=-Ext_Jacob\dM;
    
    
    alpha=alpha+rr*corr(1,1);
    V2=V2+rr*corr(2,1);
    Lambda=Lambda+rr*corr(3,1);
    hat_Lambda=hat_Lambda+rr*corr(4,1);
    
    Iter = Iter + 1;
    Tole=max(abs(dM));  
    counter = counter + 1;
    if counter ==150;
        break;
    end    
end


VV2(Iterr)=V2;

Lambda2(Iterr)=hat_Lambda*Ql;
Lambda1(Iterr)=Pl*(1+Lambda);

Tole_corr(Iterr)=Tole;


Iterr = Iterr + 1;
    Tolee=max(abs(scale));  
    counterr = counterr + 1;
    if counterr ==40;
        break;
    end    
end

LL=[-Lambda1 0 Lambda1];
VV=[VV2 0.5 VV2];

plot(LL,VV,'+')
%plot(Lambda1,VV2,'+')
xlim([-Inf Inf])
ylim([-inf Inf])


