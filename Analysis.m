function Model=Analysis(Model,Force,U_s,Analysis_Options)


disp('Analysis');
disp('-----------------------');
%----------


if Model.Analysis_Status.Analyzed==0
Model=Model_Set_Initial_Conditions(Model);                %Initial Conditions
Modeli=Model;                                             %Intitial Model       
Model=Model_Update_Elements(Model,Modeli);disp('-Update Model Elements');       
else
Modeli=Model;                                             %Intitial Model       
end

[Support]=Get_Support(Model);                             %Get Constraint Vector
[U_s_i]=Get_U_s(Model);                                   %Get Initial Constraint Displacement
if isempty(U_s);  U_s=U_s_i; end                          %Use initial Constraint Dispalcement (as applied constraint disp) if it is empty
Model=Set_U_s(Model,U_s);                                 %Update Model Constraint Displacement
dU_s=U_s-U_s_i;                                           %Incremental Constraint Displacement


NDOF=Model.NDOF;
Free=[1:1:NDOF]'; Free(Support)=[];                       %Get Free DOF
DU=zeros(NDOF,1); DU(Support)= dU_s;                      %Initial Displacement Increment


Model.Force=Force;                                        %Set Model external Force
F_ex=Get_F_ex(Force,NDOF);                                %Get external force Vector


U=Get_U(Model);                                           %Get initial displacement from model 

Model.Analysis_Status.Convergance=0;

%Start Iterations
%------------------------------------

iter=0;                                                          %Initial Iteration number
[F_i]=Model_Get_F_i(Model);disp('-Get Model Internal Forces');   %Get Initial internal force from model 
dF=F_ex-F_i;                                                     %Initial Residual Force
Diverge=0;                                                       %Scalar used to detect divergence in iterations


while 1


disp(['Iteration ' num2str(iter)]);


U_i=U;                                                          %Displacement from previous iteration


if strcmp(Analysis_Options.Iteration_Method,'Newton-Raphson') ||...
  (strcmp(Analysis_Options.Iteration_Method,'Modified Newton-Raphson') && iter==0)        
K=Global_Stiffness(Model);                                      %Get global stiffness
[K_ff, K_fs]=Stiffness_Submatrices(K,Support,Free);             %Matrix Decomposition
end


dF_eff=dF(Free)-K_fs*dU_s;                               %Effective Force Increment  
if iter==0; dF_eff0=dF_eff; end
[dU_f]=Solve_Linear_System(K_ff,dF_eff);                 %Solve Linear System


eta=1;                                                   %Line Search Method  
if Analysis_Options.Line_Search_Used==1   
if iter>0 || Analysis_Options.Line_Search_In_Iter_0==1
eta=Line_Search(Model,Modeli,U_i,dU_f,Free,dU_s,Support,F_ex,Analysis_Options.Line_Search_Options);
if isnan(eta); break; end
end
end
dU(Free,1)=eta*dU_f; dU(Support,1)=dU_s;                 %Displacement Incremenet


DU=DU+dU;                                                          %Total displacement increment
U=U_i+dU;                                                          %New displacement
Model=Set_U(Model,U);disp('-Update Model Displacement');           %Update Model displacement
Model=Model_Update_Elements(Model,Modeli);disp('-Update Model Elements'); %Update Stiffness after displacement update  
[F_i]=Model_Get_F_i(Model);disp('-Get Model Internal Forces');     %Get internal force from model 
dF=F_ex-F_i;                                                       %Residual force
dF_f=dF(Free);                                                     %Free DOF Residual Forces                            


e_F=norm(dF_f);                                                    %Force norm ratio 
if Analysis_Options.Tolerance_Relative==1; e_F=e_F/norm(dF_eff0); end                                   
disp(['-Force norm='  num2str(e_F)]);
%----
if iter>0
e_U=norm(dU);                                                      %Displacement norm ratio 
if Analysis_Options.Tolerance_Relative==1; e_U=e_U/norm(dU0); end                                   
disp(['-Disp norm='  num2str(e_U)]);
end
%----

e_F_iter(iter+1)=e_F;                                              %Save iteration parameters in vectors
if iter>0; e_U_iter(iter+1)=e_U; end  
Model_iter{iter+1}=Model;
F_i_iter{iter+1}=F_i;


if e_F<Analysis_Options.Force_Tolerance;  disp('***Iterations Converged due to Force Tolerance');
    Model.Analysis_Status.Convergance=1; break; end
if iter>0; if e_U<Analysis_Options.Disp_Tolerance; disp('***Iterations Converged due to Displacement Tolerance');
        Model.Analysis_Status.Convergance=1; break; end; end
if iter==Analysis_Options.Max_Iterations;  disp('***Max Number of Iterations Reached'); break; end


if iter>0 && e_F_iter(iter+1)>e_F_iter(iter); Diverge=Diverge+1; else; Diverge=0; end
if Diverge>=Analysis_Options.Diverge_Tolerance; disp('***Iterations Terminated Due to Divergance'); break; end


if iter==0 
dU_s(1:end)=0;                                                  %Incremental Support displacement vanishes after iteration zero
dU0=dU;                                                         %Initial Incremental displacement
end


iter=iter+1;                                                    %Iteration number


end

%End Iterations
%------------------------------------------------

I=iter+1;
if Model.Analysis_Status.Convergance==0
if I>1; [~,I] = min(e_F_iter); end                        %Find Iteration with minimium residual force
if I~=iter+1; F_i=F_i_iter{I}; Model=Model_iter{I}; end  
disp(['***Iteration ' num2str(I-1) ' with Minimum Residual Force was used']);
end

F_s=Get_F_s(F_ex,F_i,Support);                            %Get final constraint force
Model=Set_F_s(Model,F_s);                                 %Update Model F_s



if Model.Analysis_Status.Analyzed==0
Model.Analysis_Status.Analyzed=1;
end
if exist('e_F_iter','var'); Model.Analysis_Status.e_F=e_F_iter(I); end
if exist('e_U_iter','var'); Model.Analysis_Status.e_U=e_U_iter(I); end  
Model.Analysis_Status.Internal_Force=F_i;
Model.Analysis_Status.Residual_Force=dF_f;
Model.Analysis_Status.Displacement_Increment=DU;
Model.Analysis_Status.Analysis_Options=Analysis_Options;

disp('  ');


end

%Tips to decrease computation time: 
%i)Internal Forces (F_i) is computed at the begining of the algorithm. it could be called from the end of previous analysis instead.
%ii) Line Search Algorithm finds (F_i) at eta=0. It could be called from the main subroutine instead.
%iii) After doing Line Search (F_i) is computed again for the chosen eta value. Instead, it could be called from Line Search Subroutine.