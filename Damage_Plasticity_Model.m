function [Material_State,D]=Damage_Plasticity_Model(Material,Material_State,e)

%Notes: 
%1-Denominator-layout notation was used for matrix calculus
%2-d_f_d_s= TT'*d_f_d_S is total derivative while d_d_f_d_d_s=TT'*d_d_f_d_d_S*TT is not total because principal stress derivative with respect to direction is zero and
%second order derivative does not equal to zero. This causes the tangent stiffness matrix not to be exact. It can be solved by calculating d_d_f_d_d_s directly
%3-Sometimes when both Rankine and Drucker-Prager yield functions are activated, the solution does not converge and then only one of them converge after assigning them.

%Material Properties
%--------------------
E=Material.E;
v=Material.v;
f_t=Material.f_t;
g_f=Material.g_f;
f_c=Material.f_c;
f_c2=Material.f_c2;

beta=sqrt(3)*(f_c2-f_c)/(2*f_c2-f_c);
Hp=f_c2*f_c/(sqrt(3)*(2*f_c2-f_c));

%Values from Previous Increment
%-------------------------------
e_i=Material_State.e;
s_i=Material_State.s_eff;   %Effective stress
k_RK_i=Material_State.k_RK;
k_DP_i=Material_State.k_DP;
k_D=Material_State.k_D;     %Doesn't require k_D_i because it is not in return mapping algorithm iterations  

%Constitutive law of elastic material
%----------------------------------------

E_hat=E/((1-2*v)*(1+v));
G=1/2*E/(1+v);

D_e= [E_hat*(1-v) E_hat*v     E_hat*v     0 0 0;...
      E_hat*v     E_hat*(1-v) E_hat*v     0 0 0;... 
      E_hat*v     E_hat*v     E_hat*(1-v) 0 0 0;... 
      0           0           0           G 0 0;...
      0           0           0           0 G 0;...
      0           0           0           0 0 G];

%Strain Increment and Elastic Stress
%-----------------------------------
d_e=e-e_i;
s_e=s_i+D_e*d_e;

%Initial Solution
%------------------
i=0;
s=s_e;
d_k_RK=0;
d_k_DP=0;

k_RK=k_RK_i;
k_DP=k_DP_i;

while 1

Converged=0;

i=i+1;

%Principal Stress
%------------------
[T,S] = eig([s(1) s(4) s(6);...
             s(4) s(2) s(5);...
             s(6) s(5) s(3)]);
S=diag(S);                 
[S,I]=sort(S,'descend'); T=T(:,I);         

TT=[T(1,1)^2 T(2,1)^2 T(3,1)^2 2*T(1,1)*T(2,1) 2*T(2,1)*T(3,1) 2*T(3,1)*T(1,1);...             %Transformation of Stress Matrix
    T(1,2)^2 T(2,2)^2 T(3,2)^2 2*T(1,2)*T(2,2) 2*T(2,2)*T(3,2) 2*T(3,2)*T(1,2);...
    T(1,3)^2 T(2,3)^2 T(3,3)^2 2*T(1,3)*T(2,3) 2*T(2,3)*T(3,3) 2*T(3,3)*T(1,3)];

%Calculate yield function and its derivatives (Rankine)
%--------------------------------------------------------
s_h_RK=f_t+0*k_RK;                                                                             %Hardening Law (No Hardening)
h_RK=0;                                                                                        %Hardening Law Derivative
d_f_RK_d_k=-h_RK;
f_RK=0;
if S(1)>0 && S(2)<=0 && S(3)<=0
f_RK=S(1)-s_h_RK;                                                                              %Rankine Yield Surface
 
d_f_RK_d_S=[1;...                                                                              %Yield function derivative
            0;...
            0];
d_d_f_RK_d_d_S=[0 0 0;...
                0 0 0;...
                0 0 0];    
        
elseif S(1)>0 && S(2)>0 && S(3)<=0
f_RK=sqrt(S(1)^2+S(2)^2)-s_h_RK;                                                               %Rankine Yield Surface

d_f_RK_d_S=[S(1)/(S(1)^2 + S(2)^2)^(1/2);...                                                   %Yield function derivative
            S(2)/(S(1)^2 + S(2)^2)^(1/2);...
            0];
d_d_f_RK_d_d_S=[ 1/(S(1)^2 + S(2)^2)^(1/2) - S(1)^2/(S(1)^2 + S(2)^2)^(3/2),  -(S(1)*S(2))/(S(1)^2 + S(2)^2)^(3/2),0;...
                  -(S(1)*S(2))/(S(1)^2 + S(2)^2)^(3/2), 1/(S(1)^2 + S(2)^2)^(1/2) - S(2)^2/(S(1)^2 + S(2)^2)^(3/2),0;...
                0 0 0];               

elseif S(1)>0 && S(2)>0 && S(3)>0 
f_RK=sqrt(S(1)^2+S(2)^2+S(3)^2)-s_h_RK;                                                        %Rankine Yield Surface 

d_f_RK_d_S=[S(1)/(S(1)^2 + S(2)^2+ S(3)^2)^(1/2);...                                           %Yield function derivative
            S(2)/(S(1)^2 + S(2)^2+ S(3)^2)^(1/2);...
            S(3)/(S(1)^2 + S(2)^2+ S(3)^2)^(1/2)]; 
d_d_f_RK_d_d_S=[ 1/(S(1)^2 + S(2)^2 + S(3)^2)^(1/2) - S(1)^2/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2),            -(S(1)*S(2))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2),                                 -(S(1)*S(3))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2);...
                            -(S(1)*S(2))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2), 1/(S(1)^2 + S(2)^2 + S(3)^2)^(1/2) - S(2)^2/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2),                                 -(S(2)*S(3))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2);...
                            -(S(1)*S(3))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2),                                -(S(2)*S(3))/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2), 1/(S(1)^2 + S(2)^2 + S(3)^2)^(1/2) - S(3)^2/(S(1)^2 + S(2)^2 + S(3)^2)^(3/2)];      

end

%Calculate yield function and its derivatives (Drucker-Prager)
%--------------------------------------------------------------
s_h_DP=Hp+0*k_DP;                                                                                        %Hardening Law (No Hardening)
h_DP=0;                                                                                                  %Hardening Law Derivative  
d_f_DP_d_k=-h_DP;
f_DP=beta/3*(S(1)+S(2)+S(3))+1/sqrt(3)*sqrt(S(1)^2+S(2)^2+S(3)^2-S(1)*S(2)-S(2)*S(3)-S(3)*S(1))-s_h_DP;  %Drucker-Prager Yield Surface

d_f_DP_d_S=[beta/3 - (3^(1/2)*(S(2) - 2*S(1) + S(3)))/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2));...%Yield function derivative
            beta/3 - (3^(1/2)*(S(1) - 2*S(2) + S(3)))/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2));...
            beta/3 - (3^(1/2)*(S(1) + S(2) - 2*S(3)))/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2))];
d_d_f_DP_d_d_S=[   3^(1/2)/(3*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(2) - 2*S(1) + S(3))^2)/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)), - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) - 2*S(2) + S(3))*(S(2) - 2*S(1) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)), - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) + S(2) - 2*S(3))*(S(2) - 2*S(1) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2));...
           - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) - 2*S(2) + S(3))*(S(2) - 2*S(1) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)),                  3^(1/2)/(3*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) - 2*S(2) + S(3))^2)/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)), - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) + S(2) - 2*S(3))*(S(1) - 2*S(2) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2));...
           - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) + S(2) - 2*S(3))*(S(2) - 2*S(1) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)), - 3^(1/2)/(6*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) + S(2) - 2*S(3))*(S(1) - 2*S(2) + S(3)))/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2)),                  3^(1/2)/(3*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(1/2)) - (3^(1/2)*(S(1) + S(2) - 2*S(3))^2)/(12*(S(1)^2 - S(1)*S(2) - S(1)*S(3) + S(2)^2 - S(2)*S(3) + S(3)^2)^(3/2))];      
                                                                      
if i==1; Yieldfn=[]; end   
if i==1 && f_RK<=0 && f_DP<=0; break; end
if i==1001 || i==2001 || i==3001 || i==4001 || i==5001 || i==6001  %Restart all the iterations and try new yield functions (for i>4000 the yield functions are fixed)
s=s_e;
d_k_RK=0;
d_k_DP=0;
k_RK=k_RK_i;
k_DP=k_DP_i;
if i==1001 || i==4001; Yieldfn=[]; Yieldfn{1}='RK'; continue; end
if i==2001 || i==5001; Yieldfn=[]; Yieldfn{1}='DP'; continue; end
if i==3001 || i==6001; Yieldfn=[]; Yieldfn{1}='RK'; Yieldfn{2}='DP'; continue; end
end

%Set of Initial Active Yield Functions
%--------------------------------------
if i==1   
c=0;
if f_RK>0 
c=c+1;    
Yieldfn_i{c}='RK';
end
if f_DP>0 
c=c+1;    
Yieldfn_i{c}='DP';
end
Yieldfn=Yieldfn_i;    %Set of Active Yield Functions 
end

while 1
%Use only Active Yield Functions
%--------------------------------
c=0;
f=[]; h=[]; d_f_d_s=[]; d_d_f_d_d_s=[]; d_f_d_k=[]; k_i=[]; d_k=[];
if not(isempty(find(contains(Yieldfn,'RK'), 1))) 
c=c+1;
f(c,1)=f_RK;
h(c,c)=h_RK;
d_f_RK_d_s=TT'*d_f_RK_d_S;
d_d_f_RK_d_d_s=TT'*d_d_f_RK_d_d_S*TT;
d_f_d_s(:,c)=d_f_RK_d_s;
d_d_f_d_d_s{c}=d_d_f_RK_d_d_s;
d_f_d_k(c,c)=d_f_RK_d_k;
k_i(c,1)=k_RK_i;
d_k(c,1)=d_k_RK;
end
if not(isempty(find(contains(Yieldfn,'DP'), 1)))   
c=c+1;
f(c,1)=f_DP;
h(c,c)=h_DP;
d_f_DP_d_s=TT'*d_f_DP_d_S;
d_d_f_DP_d_d_s=TT'*d_d_f_DP_d_d_S*TT;
d_f_d_s(:,c)=d_f_DP_d_s;
d_d_f_d_d_s{c}=d_d_f_DP_d_d_s;
d_f_d_k(c,c)=d_f_DP_d_k;
k_i(c,1)=k_DP_i;
d_k(c,1)=d_k_DP;
end              

%Calculate Risidue
%-------------------
r1=D_e\(s-s_e);
for j=1:1:c
r1=r1+d_k(j)*d_f_d_s(:,j);    
end
r2=f;
r=[r1 ; r2];

if norm(r1)<10^(-6)*f_t && norm(r2)<10^(-6)*f_t; Converged=1; break; end  

%Calculate Jacobian
%--------------------
J11=inv(D_e);
for j=1:1:c
J11=J11+d_k(j)*d_d_f_d_d_s{j};    
end
J12=d_f_d_s;   
J21=d_f_d_s';  
J22=d_f_d_k;
J=[J11 J12 ; J21 J22];

%Calculate Variation of Stress Vector and Hardening Parameter
%--------------------------------------------------------------
dd_sk=-J\r; 
dd_s=dd_sk(1:6);
dd_k=dd_sk(7:end);
if i<4000 && not(isempty(find(d_k+dd_k<0, 1))); Yieldfn(d_k+dd_k<0)=[]; continue; end   %Check if there is negative dk, remove the yield function and restart the iteration 
s=s+dd_s;
d_k=d_k+dd_k;
k=k_i+d_k;

for j=1:1:c
if Yieldfn{j}=='RK'; d_k_RK=d_k(j); k_RK=k(j); end
if Yieldfn{j}=='DP'; d_k_DP=d_k(j); k_DP=k(j); end
end

break;

end

if Converged==1; break; end

end

%Calculate Tangential Elastoplastic Stiffness Matrix
%-----------------------------------------------------

if isempty(Yieldfn)              %There is no active yield function: D_ep=D_e 

D_ep=D_e;
d=1-exp(-f_t*k_D/g_f);
s_eff=s;                         %Stress from return mapping algorithm is the effective stress
s=(1-d)*s_eff;
D=(1-d)*D_ep;

else
    
Hinv=inv(D_e);
for j=1:1:c
Hinv=Hinv+d_k(j)*d_d_f_d_d_s{j};    
end    
H=inv(Hinv);
D_ep=H-H*d_f_d_s*inv(h+d_f_d_s'*H*d_f_d_s)*d_f_d_s'*H;

%Plastic Strain and its Derivative
%----------------------------------
e_p=e-D_e\s;
d_e_p_d_e=inv(D_e)*(D_e-D_ep);

%Principal Plastic Strain
%------------------------------
[Tp,Ep] = eig([e_p(1) e_p(4)/2 e_p(6)/2;...
               e_p(4)/2 e_p(2) e_p(5)/2;...
               e_p(6)/2 e_p(5)/2 e_p(3)]);
Ep=diag(Ep);                 
[Ep,I]=sort(Ep,'descend'); Tp=Tp(:,I);         

TTp=[Tp(1,1)^2 Tp(2,1)^2 Tp(3,1)^2 Tp(1,1)*Tp(2,1) Tp(2,1)*Tp(3,1) Tp(3,1)*Tp(1,1);...  %Transformation of Plastic Strain Matrix
     Tp(1,2)^2 Tp(2,2)^2 Tp(3,2)^2 Tp(1,2)*Tp(2,2) Tp(2,2)*Tp(3,2) Tp(3,2)*Tp(1,2);...
     Tp(1,3)^2 Tp(2,3)^2 Tp(3,3)^2 Tp(1,3)*Tp(2,3) Tp(2,3)*Tp(3,3) Tp(3,3)*Tp(1,3)];

%Equivelent Plastic Strain and Its Derivative 
%---------------------------------------------- 
Ep_eqv=sqrt(Ep(1)^2+Ep(2)^2+Ep(3)^2);

d_Ep_eqv_d_Ep=[Ep(1)/(Ep(1)^2 + Ep(2)^2+ Ep(3)^2)^(1/2);...                                           
               Ep(2)/(Ep(1)^2 + Ep(2)^2+ Ep(3)^2)^(1/2);...
               Ep(3)/(Ep(1)^2 + Ep(2)^2+ Ep(3)^2)^(1/2)];

d_Ep_eqv_d_ep=TTp'*d_Ep_eqv_d_Ep;

d_Ep_eqv_d_e=d_e_p_d_e'*d_Ep_eqv_d_ep;

%Loading and Unloading Conditions
%---------------------------------
if Ep_eqv>k_D
k_D=Ep_eqv;    
d_k_D_d_Ep_eqv=1;
else
d_k_D_d_Ep_eqv=0;    
end

d_k_D_d_e=d_k_D_d_Ep_eqv*d_Ep_eqv_d_e;

%Damage and its derivative
%--------------------------
d=1-exp(-f_t*k_D/g_f);

d_d_d_k_D=f_t/g_f*exp(-f_t*k_D/g_f);

d_d_d_e=d_d_d_k_D*d_k_D_d_e;

%Stress and its derivative
%---------------------------
s_eff=s;         %Stress from return mapping algorithm is the effective stress
s=(1-d)*s_eff;
D=(1-d)*D_ep-s_eff*d_d_d_e';

end

%Final Material State
%---------------------
Material_State.e=e;
Material_State.s=s;
Material_State.s_eff=s_eff;
Material_State.k_RK=k_RK;
Material_State.k_DP=k_DP;
Material_State.k_D=k_D;
Material_State.d=d;

end

