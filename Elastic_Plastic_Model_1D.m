function [Material_State,D]=Elastic_Plastic_Model_1D(Material,Material_State,e)

%Material properties
%--------------------
E=Material.E;
f_y=Material.f_y;

%Values from previous step
%-------------------------
e_i=Material_State.e;
s_i=Material_State.s;
k_i=Material_State.k;

%Strain increment and elastic stress
%-----------------------------------
d_e=e-e_i;
s_e=s_i+E*d_e;    

%Calculate yield function 
%-------------------------
s_h=f_y+0*E*k_i;  %Hardening Law 
h=0*E; %Hardening Law Derivative
f=abs(s_e)-s_h;

if f<0
s=s_e;
k=k_i;
D=E;
else
%Return mapping
%----------------
dk=f/(E+h);
s=(1-dk*E/abs(s_e))*s_e;
k=k_i+dk;
D=E*h/(E+h);
end

Material_State.e=e;
Material_State.s=s;
Material_State.k=k;

end
