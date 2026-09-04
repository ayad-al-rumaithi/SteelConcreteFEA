function [Material_State,D]=Traction_Seperation_Model_1D(Material,Material_State,e)

%Material properties
%--------------------
Tao_0=Material.Tao_0;
delta_0=Material.delta_0;
delta_max=Material.delta_max;

%Values from previous step
%-------------------------
k_i=Material_State.k;

Kt=Tao_0/delta_0;
Ks=-Tao_0/(delta_max-delta_0);
e_abs=abs(e);
   
if e_abs>=k_i   %Loading
[s, D]=bilinear(delta_0,delta_max,Kt,Ks,e);
k=e_abs;
else            %Unloading/Reloading
[smax, ~]=bilinear(delta_0,delta_max,Kt,Ks,k_i);
s=abs(smax)/k_i*e;
D=abs(smax)/k_i; 
k=k_i;
end
    
Material_State.e=e;
Material_State.s=s;
Material_State.k=k;

end

function [s, D]=bilinear(delta_0,delta_max,Kt,Ks,e)

e_abs=abs(e);
e_sign=sign(e);
    
if e_abs <=delta_0
s=e_sign*Kt*e_abs;
D=Kt;    
elseif e_abs <=delta_max
s=e_sign*(Kt*delta_0+Ks*(e_abs-delta_0));    
D=Ks;    
else
s=0;    
D=0;    
end 

end
