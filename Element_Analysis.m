function [Element]=Element_Analysis(Element,Material,V,U)

Element_Type=Element.Type;
Integration_Points=Element.Integration_Points;
UU=reshape(U',size(U,1)*size(U,2),1);
NUU=length(UU);
N_Int=Element.No_Int_Points;
J_det=zeros(N_Int,1);
B=cell(N_Int,1);
G=cell(N_Int,1);
K=zeros(NUU,NUU);
f=zeros(NUU,1);
Material_hat=Material;
%--------------------------------------------------------------------------
if strcmp(Element_Type,'wdg6')              %Wedge Element
Element_Category='Solid';    
F_eq=2;    
if N_Int==6
r=[0.166666667;  0.666666667;  0.166666667;  0.166666667; 0.666666667; 0.166666667];
s=[0.166666667;  0.166666667;  0.666666667;  0.166666667; 0.166666667; 0.666666667];
t=[-0.577350269; -0.577350269; -0.577350269; 0.577350269; 0.577350269; 0.577350269];
w=[0.166666667;  0.166666667;  0.166666667;  0.166666667; 0.166666667; 0.166666667];
end
for i=1:1:N_Int
[J_det(i), B{i}]=wdg6(V,r(i),s(i),t(i));
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Type,'tetra4')            %Tetrahedral Element
Element_Category='Solid';    
F_eq=6;    
if N_Int==1
r=[1/4];
s=[1/4];
t=[1/4];
w=[1/6];
end
for i=1:1:N_Int
[J_det(i), B{i}]=tetra4(V,r(i),s(i),t(i));
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Type,'brick8')            %Brick Element
Element_Category='Solid';    
F_eq=1;    
if N_Int==8
r=1/sqrt(3)*[-1; 1; 1; -1; -1; 1; 1; -1];
s=1/sqrt(3)*[-1; -1; 1; 1; -1; -1; 1; 1];
t=1/sqrt(3)*[-1; -1; -1; -1; 1; 1; 1; 1];
w=[1; 1; 1; 1; 1; 1; 1; 1];
end
for i=1:1:N_Int
[J_det(i), B{i}]=brick8(V,r(i),s(i),t(i));
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Type,'brick8 TL')         %Brick Element (Total Lagrangian)
Element_Category='Solid TL';    
F_eq=1;    
if N_Int==8
r=1/sqrt(3)*[-1; 1; 1; -1; -1; 1; 1; -1];
s=1/sqrt(3)*[-1; -1; 1; 1; -1; -1; 1; 1];
t=1/sqrt(3)*[-1; -1; -1; -1; 1; 1; 1; 1];
w=[1; 1; 1; 1; 1; 1; 1; 1];
end
for i=1:1:N_Int
[J_det(i), B{i}, G{i}]=brick8_TL(V,UU,r(i),s(i),t(i));
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Type,'truss2')            %Truss Element
Element_Category='Bar';    
F_eq=1;    
if N_Int==1
r=[0];
w=[2];
end
for i=1:1:N_Int
[J_det(i), B{i}]=truss2(V,r(i));
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Type,'spring2')           %Spring Element
Element_Category='Spring';    
if N_Int==1
r=[0];
w=[1];
end
Direction=Element.Direction; 
for i=1:1:N_Int
[J_det(i), B{i}]=spring2(V,Direction,r(i));
end
end
%--------------------------------------------------------------------------
for i=1:1:N_Int                             %Check Jacobain
if J_det(i)<0; error('Negative Jacobian'); end
end
%--------------------------------------------------------------------------
if strcmp(Element_Category,'Solid')...      %Solid Elements
|| strcmp(Element_Category,'Solid TL')       
Vol=w'*J_det;
Element.Vol=Vol;
if isfield(Material,'G_f')   %use normalized fracture energy
l_eq=(F_eq*Vol)^(1/3);
Material_hat=rmfield(Material_hat,'G_f');
Material_hat.g_f=Material.G_f/l_eq;
end
for i=1:1:N_Int   
e=B{i}*UU; 
[Integration_Points{i},D]=Material_Model(Material_hat,Integration_Points{i},e);
K=K+w(i)*J_det(i)*B{i}'*D*B{i};
s=Integration_Points{i}.s;
f=f+w(i)*J_det(i)*B{i}'*s;
Integration_Points{i}.B=B{i};
Integration_Points{i}.J_det=J_det(i);
Integration_Points{i}.D=D;
end
end
%--------------------------------------------------------------------------
if strcmp(Element_Category,'Solid TL')      %Solid Elements (Total Lagrangian)
z1=zeros(3,3);    
for i=1:1:N_Int 
s=Integration_Points{i}.s;
s1=[s(1) s(4) s(6);...
    s(4) s(2) s(5);...
    s(6) s(5) s(3)];
S=[s1 z1 z1;...
   z1 s1 z1;... 
   z1 z1 s1];      
K=K+w(i)*J_det(i)*G{i}'*S*G{i};
Integration_Points{i}.G=G{i};
end    
end
%--------------------------------------------------------------------------
if strcmp(Element_Category,'Bar')           %Bar Elements
Length=w'*J_det;
Element.Length=Length;
if isfield(Material,'G_f')   %use normalized fracture energy
l_eq=F_eq*Length;
Material_hat=rmfield(Material_hat,'G_f');
Material_hat.g_f=Material.G_f/l_eq;
end
A=Element.A; %section area
for i=1:1:N_Int   
e=B{i}*UU; 
[Integration_Points{i},D]=Material_Model(Material_hat,Integration_Points{i},e);
K=K+A*w(i)*J_det(i)*B{i}'*D*B{i};
s=Integration_Points{i}.s;
f=f+A*w(i)*J_det(i)*B{i}'*s;
Integration_Points{i}.B=B{i};
Integration_Points{i}.J_det=J_det(i);
Integration_Points{i}.D=D;
end    
end
%--------------------------------------------------------------------------
if strcmp(Element_Category,'Spring')        %Spring Elements
for i=1:1:N_Int   
e=B{i}*UU; 
[Integration_Points{i},D]=Material_Model(Material_hat,Integration_Points{i},e);
K=K+w(i)*J_det(i)*B{i}'*D*B{i};
s=Integration_Points{i}.s;
f=f+w(i)*J_det(i)*B{i}'*s;
Integration_Points{i}.B=B{i};
Integration_Points{i}.J_det=J_det(i);
Integration_Points{i}.D=D;
end 
end
%--------------------------------------------------------------------------
Element.K=K;
Element.f=f;
Element.Integration_Points=Integration_Points;

end
