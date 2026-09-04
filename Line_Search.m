function eta=Line_Search(Model,Modeli,U_i,dU_f,Free,dU_s,Support,F_ex,Line_Search_Options)

%This algorithm was taken from Section 9.2 in
%Crisfield M.A. Vol.1. Non-Linear Finite Element Analysis of Solids and Structures.. Essentials 
%and then modified by user (Added Negative Line Search)
%----------------------

sign=1; %Posative line search

disp('-Line Search');

x_min=Line_Search_Options.x_min;
x_max=Line_Search_Options.x_max;
beta=Line_Search_Options.beta;
Max_Iter=Line_Search_Options.Max_Iter;

amp=2;

x(1)=0; s(1)=Get_S(Model,Modeli,U_i,dU_f,Free,dU_s,Support,F_ex,x(1)); r(1)=1;

x(2)=1; 

i=1;

while 1

s(i+1)=Get_S(Model,Modeli,U_i,dU_f,Free,dU_s,Support,F_ex,x(i+1));
    
r(i+1)=s(i+1)/s(1);

disp(['--eta=' num2str(x(i+1)) ' r=' num2str(r(i+1))]);

if abs(r(i+1))<beta; break; end                                              %Convergence criteria

if i==Max_Iter; disp('--***Max Number of Iterations Reached'); break; end    %Max number of iterations

if x(i+1)<0 && sign==1 
sign=-1; 
disp('--***Switch to Negative Line Search'); 
x=x(1); r=r(1); i=1; x(2)=-1;
continue;
end
    
x(i+2)=Search(x,r,x_min,x_max,amp,sign);

i=i+1;

if x(i+1)==x(i); disp('--***Iterations Terminated'); x(i+1)=[]; break; end   %Convergence stopped

end

[~,I] = min(abs(r(2:end)));

eta=x(I+1);

disp(['--Optimum eta=' num2str(eta)]);

end


function x_new=Search(x,r,x_min,x_max,amp,sign)

if isempty(find(r<0, 1))

x_new=x(end)-(x(end)-x(end-1))/(r(end)-r(end-1))*r(end);

if sign==1
x_max_p=max(x); if  x_new>amp*x_max_p; x_new=amp*x_max_p; end 
else
x_min_p=min(x); if  x_new<amp*x_min_p || x_new>0; x_new=amp*x_min_p; end 
end  %Avoid dangerous extrapolation 

else
    
x_l=Inf; x_r=-Inf; 

for i=1:1:length(r)
 
if x(i)<x_l && r(i)*sign<0; x_l=x(i);  r_l=r(i);  end
    
if x(i)>x_r && r(i)*sign>0; x_r=x(i);  r_r=r(i);  end       
   
end

x_new_1=x_l-(x_l-x_r)/(r_l-r_r)*r_l;

x_new_2=x_r+(0.5-sign*0.3)*(x_l-x_r);                                 %give factor 0.2 for posative and 0.8 for negative

if abs(x_new_1)>abs(x_new_2); x_new=x_new_1; else; x_new=x_new_2; end %Choose max absolute

end

if x_new>x_max;  x_new=x_max; end                                     %Check with bounadries
if x_new<-x_max; x_new=-x_max; end

if x_new<=x_min && x_new>=0;  x_new=x_min; end
if x_new>=-x_min && x_new<0;  x_new=-x_min; end

    
end


function s=Get_S(Model,Modeli,U_i,dU_f,Free,dU_s,Support,F_ex,eta)

dU(Free,1)=eta*dU_f; dU(Support,1)=dU_s;  
U=U_i+dU; 

Model1=Set_U(Model,U);  
Model1=Model_Update_Elements(Model1,Modeli); 

[F_i]=Model_Get_F_i(Model1);  
dF=F_ex-F_i; 
dF_f=dF(Free); 

s=dF_f'*dU_f;

end