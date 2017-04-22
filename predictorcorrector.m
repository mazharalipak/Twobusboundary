clc;
clear all;

%% Input Data.....
tic;
rr=1;              % alpha.....
tau=-0.05;           % scaling for continuation.......
scale=1;           % scaling for Ql.........

%% Input data got from first freeing only Pl as free parameter.....

V2=0.5;         % V2........
alpha=6.8661e-04;      % Alpha........
Lambda=-1.0034;   % Lambda 1.....

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
Ql=scale*(-0.25);

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


zp=[alpha;V2;Lambda;hat_Lambda];
%zp=z1+tau*tang_vector;                      % Predicted step on the boundry ..........


%% Corrector...................

Tole = 1;  
Iter = 1;
counter = 0;
while (Tole > 1e-6)
    
    ff1=((V2*sin(alpha)))+Pl*(1+Ldfactor*Lambda);
    ff2=(((V2^2) -V2*cos(alpha)))- hat_Lambda*Ql*(1);
    ff3=2*V2^2 *cos(alpha) -V2;
    ff4=(zp-zp)'*tang_vector-0.0;
    
    Ext_Jacob=[V2*cos(alpha) sin(alpha) Pl*Ldfactor 0; V2*sin(alpha)  2*V2-cos(alpha) 0 -Ql; -2*(V2^2)*sin(alpha) (4*V2*cos(alpha))-1 0 0;tang_vector'];
   
    
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
    if counterr ==80;
        break;
    end    
end

LL1=[0 -Lambda1];
VV1=[0.5 VV2];

LL2=[0 Lambda1];
VV3=[0.5 VV2];

%% Analytical solution of Outer Boundary of Solution Space......

Q_f=-0.25:0.1:0.7;
P_f=sqrt(1+4*Q_f)/2;
V_f=sqrt((2*Q_f+1)/2);

plot(LL1,VV1,LL2,VV3,P_f,V_f,'*',-P_f,V_f,'+')
title('Solution space P-V plane, Q as free parameter');
xlabel('P (MW, pu)');
ylabel('V (pu)');
legend('','','Analytical','Analytical')
ylim([-Inf Inf])
xlim([-Inf Inf])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])

print -djpeg filename.jpg -r600

toc;

