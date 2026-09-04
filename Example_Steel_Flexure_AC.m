clear all; clc; close all;

%Input
%------------
B=200;   %Width
H=200;   %Height
nx=5; 
nz=5;
L=2000;  %Length
t=100;

y=[0 L];

dx=B/nx;
dz=H/nz;

[xx, yy]=meshgrid([-B/2:dx:B/2],[-H/2:dz:H/2]);
xxx=xx(:);
yyy=yy(:);
Nodes=[xxx yyy];

Triangles=delaunay(Nodes(:,1),Nodes(:,2));
[Triangles, I1]=sort(Triangles,2);
[junk, I2]=sort(Triangles(:,1));
Triangles=Triangles(I2,:);
I1=I1(I2,:);

for i=1:1:size(Triangles,1)
Triangles(i,:)=Triangles(i,I1(i,:));    
end

[Nodes3D,Prisms] = Mesh2D_to_Mesh3D(Nodes,Triangles,y(1):t:y(2));
junk=Nodes3D(:,2); Nodes3D(:,2)=Nodes3D(:,3); Nodes3D(:,3)=junk; Nodes3D=Nodes3D(end:-1:1,:); %swap y with z

Material{1,1}.Name='Steel';
Material{1,1}.E=200000;                                        %Modulus of elasticity of steel
Material{1,1}.v=0.3;                                           %Poisson ratio of steel
Material{1,1}.f_y=400;                                         %Yield Stress of steel
Material{1,1}.Type='Elastic Plastic Model';

Analysis_Options.Iteration_Method='Newton-Raphson';            %Iteration Method
Analysis_Options.Max_Iterations=50;                            %Solver max number of iterations
Analysis_Options.Force_Tolerance=100;                          %Force Tolerance 
Analysis_Options.Disp_Tolerance=0;                             %Displacement Tolerance
Analysis_Options.Diverge_Tolerance=20;                         %Divergence Tolerance
Analysis_Options.Line_Search_Used=0;                           %Is Line Search Used
Analysis_Options.Line_Search_In_Iter_0=0;                      %Is Line Search Used in Predictor Iteration 
Analysis_Options.Tolerance_Relative=0;                         %Is Tolerance Relative or Absolute
Analysis_Options.psi=0;                                    

Analysis_Options.Line_Search_Options.x_min=0.1;                %Line Search Minimium eta               
Analysis_Options.Line_Search_Options.x_max=10;                 %Line Search Maximium eta
Analysis_Options.Line_Search_Options.beta=0.8;                 %Line Search Tolerance
Analysis_Options.Line_Search_Options.Max_Iter=10;              %Line Seach Max Number of Iterations 

Max_Step_Divisions=3;                                          %Max Number of step divisions

Show_Deformed=1;                                               %Show deformed shape
Mag_Factor=1;                                                  %Magnification factor
Result_type=13;                                                %Result type Null,ez,ey,ez,2exy,2eyz,2ezx,sz,sy,sz,sxy,syz,szx,Inelastic Parameters (0-17)

%Right hinge Constraints
[DOFx1, DOFy1, DOFz1, Node1]=find_DOF_from_Nodes(Nodes3D,[],y(2),-H/2);
M1=length(DOFy1);
for i=1:1:M1
Constraint{3*i-2,1}.DOF=DOFx1(i);     Constraint{3*i-2,1}.U=0;
Constraint{3*i-1,1}.DOF=DOFy1(i);     Constraint{3*i-1,1}.U=0;
Constraint{3*i,1}.DOF=DOFz1(i);       Constraint{3*i,1}.U=0;
end

%Left roller Constraints
[DOFx2, DOFy2, DOFz2, Node2]=find_DOF_from_Nodes(Nodes3D,[],y(1),-H/2);
M2=length(DOFy2);
for i=1:1:M2
Constraint{3*M1+2*i-1,1}.DOF=DOFx2(i);   Constraint{3*M1+2*i-1,1}.U=0;
Constraint{3*M1+2*i,1}.DOF=DOFz2(i);     Constraint{3*M1+2*i,1}.U=0;
end

%Midspan Load
[DOFx3, DOFy3, DOFz3, Node3]=find_DOF_from_Nodes(Nodes3D,[],(y(1)+y(2))/2,H/2);

%Prepare Model variable for the solver
%--------------------------------------

for i=1:1:size(Nodes3D,1)
    
Model.Node{i,1}.x=Nodes3D(i,1);
Model.Node{i,1}.y=Nodes3D(i,2);
Model.Node{i,1}.z=Nodes3D(i,3);

Model.Node{i,1}.DOFx=3*i-2;
Model.Node{i,1}.DOFy=3*i-1;
Model.Node{i,1}.DOFz=3*i;

end

for i=1:1:size(Prisms,1)
    
Model.Element{i,1}.Nodes=Prisms(i,:);

Model.Element{i,1}.Type='wdg6';
Model.Element{i,1}.No_Int_Points=6;
Model.Element{i,1}.Material=1;

end

Model.Material=Material;
Model.NDOF=3*size(Model.Node,1);
Model.Constraint=Constraint;

Model.Analysis_Status.Analyzed=0;


%External Load
%----------------------

U_s=zeros(3*M1+2*M2,1);

Load=-15*10^5;

LL=length(DOFz3);    
for j=1:1:LL  
Force_hat{j}.DOF=DOFz3(j); Force_hat{j}.Magnitude=Load/(LL-1); 
end
Force_hat{1}.Magnitude=Force_hat{1}.Magnitude/2; Force_hat{LL}.Magnitude=Force_hat{LL}.Magnitude/2;      

Arc_Length=100;

for i=1:1:5
    
Loading_Steps{i}.Force_hat=Force_hat;
Loading_Steps{i}.Arc_Length=Arc_Length;

end

%Solver
%--------------

Job=Run_Job_AC(Model,Loading_Steps,Analysis_Options,Max_Step_Divisions);

%Plot Results
%--------------

Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,8);
title('Stress in y-direction');

Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,13);
title('Hardening Parameter');

for i=1:1:length(Job)

Load(i)=0;
for j=1:1:LL
Load(i)=Load(i)+Job{i}.Model.Force{j}.Magnitude;    
end

Displacement(i)=Job{i}.Model.Node{Node3(ceil(end/2))}.Uz;
Convergance(i)=Job{i}.Model.Analysis_Status.Convergance;
e_F(i)=Job{i}.Model.Analysis_Status.e_F;
    
end

figure;

plot([0 -Displacement],[0 -Load],'-',-Displacement(Convergance==0),-Load(Convergance==0),'*r',...
                                     -Displacement(Convergance==1),-Load(Convergance==1),'*b');
xlabel('Displacement(mm)');
ylabel('Load (N)');
