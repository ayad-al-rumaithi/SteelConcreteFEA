clear all; clc; close all;

%Input
%------------
Diameter=12;     %Diameter
L=260;           %Length   
Nr=2;            %Number of Segments in r-direction
Nt=8;            %Number of Segments in t-direction 
nz=26;           %Number of Segments in z-direction
Ex=10;           %Extra Elements in top and bottom

t=L/nz;    
z=[-Ex*t L+Ex*t];
R=sqrt(pi/4*Diameter^2/(Nt*cosd(360/(2*Nt))*sind(360/(2*Nt))));
[Nodes, Triangles, Quads]=Circle_Mesh(R,Nr,Nt);
[Nodes3D,Prisms,Bricks] = Mesh2D_to_Mesh3D_2(Nodes,Triangles,Quads,z(1):t:z(2));

Material{1,1}.Name='Steel';
Material{1,1}.E=200000;                                        %Modulus of elasticity of steel
Material{1,1}.v=0.3;                                           %Poisson ratio of steel
Material{1,1}.f_y=400;                                         %Yield Stress of steel
Material{1,1}.Ep_eqv_0=0.112;
Material{1,1}.G_f=120;                                         %Gf or gf 
Material{1,1}.Type='Elastic Plastic Ductile Damage Model';

Analysis_Options.Iteration_Method='Newton-Raphson';            %Iteration Method
Analysis_Options.Max_Iterations=20;                            %Solver max number of iterations
Analysis_Options.Force_Tolerance=100;                          %Force Tolerance 
Analysis_Options.Disp_Tolerance=0;                             %Displacement Tolerance
Analysis_Options.Diverge_Tolerance=20;                         %Divergence Tolerance
Analysis_Options.Line_Search_Used=1;                           %Is Line Search Used
Analysis_Options.Line_Search_In_Iter_0=0;                      %Is Line Search Used in Predictor Iteration
Analysis_Options.Tolerance_Relative=0;                         %Is Tolerance Relative or Absolute

Analysis_Options.Line_Search_Options.x_min=0.001;              %Line Search Minimium eta               
Analysis_Options.Line_Search_Options.x_max=10;                 %Line Search Maximium eta
Analysis_Options.Line_Search_Options.beta=0.8;                 %Line Search Tolerance
Analysis_Options.Line_Search_Options.Max_Iter=10;              %Line Seach Max Number of Iterations  

Max_Step_Divisions=3;                                          %Max Number of step divisions

Show_Deformed=1;                                               %Show deformed shape
Mag_Factor=1;                                                  %Magnification factor
Result_type=13;                                                %Result type Null,ez,ey,ez,2exy,2eyz,2ezx,sz,sy,sz,sxy,syz,szx,Inelastic Parameters (0-15)

%Top vertical Constraints
for i=1:1:Nt
[DOFx, DOFy, DOFz, Node]=find_DOF_from_Nodes(Nodes3D,R*cosd(360*(i-1)/Nt),R*sind(360*(i-1)/Nt),L+(Ex-2)*t);    
Constraint{i,1}.DOF=DOFz;     Constraint{i,1}.U=0;
XX(i)=Nodes3D(Node,1);
end

%Bottom vertical Constraints
for i=1:1:Nt
[DOFx, DOFy, DOFz, Node]=find_DOF_from_Nodes(Nodes3D,R*cosd(360*(i-1)/Nt),R*sind(360*(i-1)/Nt),(-Ex+2)*t);    
Constraint{Nt+i,1}.DOF=DOFz;     Constraint{Nt+i,1}.U=0;
end

%Bottom Horizontal Constraints
[DOFx1, DOFy1, DOFz1, Node1]=find_DOF_from_Nodes(Nodes3D,R,0,-Ex*t); 
Constraint{2*Nt+1,1}.DOF=DOFx1;   Constraint{2*Nt+1,1}.U=0;
Constraint{2*Nt+2,1}.DOF=DOFy1;   Constraint{2*Nt+2,1}.U=0;

%Bottom Horizontal Constraints
[DOFx2, DOFy2, DOFz2, Node2]=find_DOF_from_Nodes(Nodes3D,-R,0,-Ex*t); 
Constraint{2*Nt+3,1}.DOF=DOFy2;   Constraint{2*Nt+3,1}.U=0;

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

NPrisms=size(Prisms,1);
for i=1:1:NPrisms
    
Model.Element{i,1}.Nodes=Prisms(i,:);
Model.Element{i,1}.Type='wdg6';
Model.Element{i,1}.No_Int_Points=6;
Model.Element{i,1}.Material=1;

end

for i=1:1:size(Bricks,1)
    
Model.Element{NPrisms+i,1}.Nodes=Bricks(i,:);
Model.Element{NPrisms+i,1}.Type='brick8';
Model.Element{NPrisms+i,1}.No_Int_Points=8;
Model.Element{NPrisms+i,1}.Material=1;

end

Model.Material=Material;
Model.NDOF=3*size(Model.Node,1);
Model.Constraint=Constraint;

Model.Analysis_Status.Analyzed=0;

%External Load
%----------------------

U_s=zeros(2*Nt+3,1);

Displacement=[1:1:4 5:5:100 102:2:130]*(L+Ex*t)/1000; 

for i=1:1:length(Displacement)

for j=1:1:Nt    
U_s(j)=Displacement(i);    
end

Loading_Steps{i}.Force=[];
Loading_Steps{i}.U_s=U_s;

end

%Solver
%--------------

Job=Run_Job(Model,Loading_Steps,Analysis_Options,Max_Step_Divisions);

%Plot Results
%--------------
Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,9);
title('Stress in z-direction');

Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,13);
title('Hardening Parameter');

Plot_Results(Job{end}.Model,Show_Deformed,Mag_Factor,15);
title('Damage');

%Deflection Nodes
[DOFxg1, DOFyg1, DOFzg1, Nodeg1]=find_DOF_from_Nodes(Nodes3D,0,0,L/2-100);
[DOFxg2, DOFyg2, DOFzg2, Nodeg2]=find_DOF_from_Nodes(Nodes3D,0,0,L/2+100);

for i=1:1:length(Job)
     
Load(i)=0;

for j=1:1:Nt
Load(i)=Load(i)+Job{i}.Model.Constraint{j}.F;
end

Displacement(i)=Job{i}.Model.Constraint{1}.U;
Deflection(i)=Job{i}.Model.Node{Nodeg2,1}.Uz-Job{i}.Model.Node{Nodeg1,1}.Uz;
Convergance(i)=Job{i}.Model.Analysis_Status.Convergance;
e_F(i)=Job{i}.Model.Analysis_Status.e_F;
    
end

figure;

plot([0 Displacement],[0 Load],'-',Displacement(Convergance==0),Load(Convergance==0),'*r',...
                                   Displacement(Convergance==1),Load(Convergance==1),'*b');
xlabel('Displacement (mm)');
ylabel('Load (N)');

figure;

plot([0 Deflection],[0 Load],'-',Deflection(Convergance==0),Load(Convergance==0),'*r',...
                                 Deflection(Convergance==1),Load(Convergance==1),'*b');
xlabel('Gauge Reading (mm)');
ylabel('Load (N)');
