clear all; clc; close all;

%Input
%------------
B=200;   %Width
H=200;   %Height
nx=5; 
ny=5; 
L=200;  %Length
t=40;

z=[0 L];

dx=B/nx;
dy=H/ny;

[xx, yy]=meshgrid([-B/2:dx:B/2],[-H/2:dy:H/2]);
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

[Nodes3D,Prisms] = Mesh2D_to_Mesh3D(Nodes,Triangles,z(1):t:z(2));

Material{1,1}.Name='Concrete';
Material{1,1}.E=30000;                                         %Modulus of elasticity of concrete
Material{1,1}.v=0.2;                                           %Poisson ratio of concrete
Material{1,1}.f_t=3.33;
Material{1,1}.G_f=0.124*2;
Material{1,1}.f_c=29; 
Material{1,1}.Type='Damage Plasticity Model';

Analysis_Options.Iteration_Method='Newton-Raphson';            %Iteration Method
Analysis_Options.Max_Iterations=10;                            %Solver max number of iterations
Analysis_Options.Force_Tolerance=1000;                         %Force Tolerance 
Analysis_Options.Disp_Tolerance=0;                             %Displacement Tolerance
Analysis_Options.Diverge_Tolerance=20;                         %Divergence Tolerance
Analysis_Options.Line_Search_Used=1;                           %Is Line Search Used
Analysis_Options.Line_Search_In_Iter_0=0;                      %Is Line Search Used in Predictor Iteration 
Analysis_Options.Tolerance_Relative=0;                         %Is Tolerance Relative or Absolute 

Analysis_Options.Line_Search_Options.x_min=0.1;                %Line Search Minimium eta               
Analysis_Options.Line_Search_Options.x_max=10;                 %Line Search Maximium eta
Analysis_Options.Line_Search_Options.beta=0.8;                 %Line Search Tolerance
Analysis_Options.Line_Search_Options.Max_Iter=10;              %Line Seach Max Number of Iterations  

Max_Step_Divisions=6;                                          %Max Number of step divisions

Show_Deformed=1;                                               %Show deformed shape
Mag_Factor=20;                                                 %Magnification factor
Result_type=15;                                                %Result type Null,ez,ey,ez,2exy,2eyz,2ezx,sz,sy,sz,sxy,syz,szx,Inelastic Parameters (0-15)

%Top vertical Constraints
[DOFx1, DOFy1, DOFz1, Node1]=find_DOF_from_Nodes(Nodes3D,[],[],z(2));
M1=length(DOFy1);
for i=1:1:M1
Constraint{i,1}.DOF=DOFz1(i);     Constraint{i,1}.U=0;
end

%Bottom vertical Constraints
[DOFx2, DOFy2, DOFz2, Node2]=find_DOF_from_Nodes(Nodes3D,[],[],z(1));
M2=length(DOFy2);
for i=1:1:M2
Constraint{M1+i,1}.DOF=DOFz2(i);   Constraint{M1+i,1}.U=0;
end

%Bottom Horizontal Constraints
[DOFx3, DOFy3, DOFz3, Node3]=find_DOF_from_Nodes(Nodes3D,B/2,H/2,z(1));
Constraint{M1+M2+1,1}.DOF=DOFx3;   Constraint{M1+M2+1,1}.U=0;
Constraint{M1+M2+2,1}.DOF=DOFy3;   Constraint{M1+M2+2,1}.U=0;

%Bottom Horizontal Constraints
[DOFx4, DOFy4, DOFz4, Node4]=find_DOF_from_Nodes(Nodes3D,-B/2,H/2,z(1));
Constraint{M1+M2+3,1}.DOF=DOFy4;   Constraint{M1+M2+3,1}.U=0;

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

U_s=zeros(M1+M2+3,1);

Displacement=-0.001*[5:5:500];

for i=1:1:length(Displacement)

U_s(1:1:M1)=Displacement(i);     

Loading_Steps{i}.Force=[];
Loading_Steps{i}.U_s=U_s;

end

%Solver
%--------------

Job=Run_Job(Model,Loading_Steps,Analysis_Options,Max_Step_Divisions);

%Plot Results
%--------------

Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,Result_type);
title('Damage');

for i=1:1:length(Job)
     
Load(i)=0;

for j=1:1:M1
Load(i)=Load(i)+Job{i}.Model.Constraint{j}.F;
end

Displacement(i)=Job{i}.Model.Constraint{1}.U;
Convergance(i)=Job{i}.Model.Analysis_Status.Convergance;
e_F(i)=Job{i}.Model.Analysis_Status.e_F;
    
end


figure;

plot([0 -Displacement],[0 -Load],'-',-Displacement(Convergance==0),-Load(Convergance==0),'*r',...
                                     -Displacement(Convergance==1),-Load(Convergance==1),'*b');
xlabel('Displacement(mm)');
ylabel('Load (N)');
