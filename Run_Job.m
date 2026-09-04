function Job=Run_Job(Model,Loading_Steps,Analysis_Options,Max_Step_Divisions)

disp('Run Job');
disp('-------------');    
disp(' ');    


Div=zeros(length(Loading_Steps),1);          %Number of divisons applied to each loading step
NDOF=Model.NDOF;

step=1;

while 1
    
if step>length(Loading_Steps); break; end
    
disp('***********************');    
disp(['Loading step:' num2str(step)]);
disp('***********************');    
disp(' ');

Analysis_Options1=Analysis_Options; 
Analysis_Options1.Force_Tolerance=Analysis_Options.Force_Tolerance; 
Analysis_Options1.Disp_Tolerance=Analysis_Options.Disp_Tolerance;   

if Analysis_Options.Tolerance_Relative==1
Analysis_Options1.Force_Tolerance=Analysis_Options1.Force_Tolerance*2^Div(step);
Analysis_Options1.Disp_Tolerance=Analysis_Options1.Disp_Tolerance*2^Div(step);
end

Model1=Analysis(Model,Loading_Steps{step}.Force,Loading_Steps{step}.U_s,Analysis_Options1); %Analysis

if Model1.Analysis_Status.Convergance==0 && Div(step)<Max_Step_Divisions
    
disp('********Loading step will be divided into two');
disp(' ');

[Loading_Steps, Div]=Divide(NDOF,Loading_Steps,step,Div);   %Divide Step into two

else

Model=Model1;    
Job{step}.Model=Model;

step=step+1;

end

end

end



function [Loading_Steps, Div]=Divide(NDOF,Loading_Steps,step,Div)

F1=zeros(NDOF,1);  
if step>1
for i=1:1:length(Loading_Steps{step-1}.Force)
    F1(Loading_Steps{step-1}.Force{i}.DOF)=Loading_Steps{step-1}.Force{i}.Magnitude;
end
end

F2=zeros(NDOF,1);                          
for i=1:1:length(Loading_Steps{step}.Force)
    F2(Loading_Steps{step}.Force{i}.DOF)=Loading_Steps{step}.Force{i}.Magnitude;
end

F=(F1+F2)/2;

Force_new=Get_Force(F);

%----------------------------------------------

if step>1
U_s_new=(Loading_Steps{step-1}.U_s+Loading_Steps{step}.U_s)/2;
else
U_s_new=Loading_Steps{step}.U_s/2;  
end

%--------------------------

Div(step)=Div(step)+1;
for i=length(Loading_Steps):-1:step
Loading_Steps{i+1}=Loading_Steps{i};
Div(i+1)=Div(i);
end
Loading_Steps{step}.Force=Force_new;
Loading_Steps{step}.U_s=U_s_new;

end