function Job=Run_Job_AC(Model,Loading_Steps,Analysis_Options,Max_Step_Divisions)

disp('Run Job');
disp('-------------');    
disp(' ');    


Div=zeros(length(Loading_Steps),1);          %Number of divisons applied to each loading step

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

Model1=Analysis_AC(Model,Loading_Steps{step}.Force_hat,Loading_Steps{step}.Arc_Length,Analysis_Options); %Analysis


if Model1.Analysis_Status.Convergance==0 && Div(step)<Max_Step_Divisions
    
disp('********Loading step will be divided into two');
disp(' ');

[Loading_Steps, Div]=Divide(Loading_Steps,step,Div);   %Divide Step into two

else

Model=Model1;    
Job{step}.Model=Model;

step=step+1;

end

end

end



function [Loading_Steps, Div]=Divide(Loading_Steps,step,Div)

Arc_Length_new=Loading_Steps{step}.Arc_Length/2;

Div(step)=Div(step)+1;
for i=length(Loading_Steps):-1:step
Loading_Steps{i+1}=Loading_Steps{i};
Div(i+1)=Div(i);
end
Loading_Steps{step}.Arc_Length=Arc_Length_new;
Loading_Steps{step+1}.Arc_Length=Arc_Length_new;

end