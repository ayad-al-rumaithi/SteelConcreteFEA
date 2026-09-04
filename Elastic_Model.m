function [Material_State,D]=Elastic_Model(Material,Material_State,e)

E=Material.E;
v=Material.v;

E_hat=E/((1-2*v)*(1+v));
G=1/2*E/(1+v);

D=[E_hat*(1-v) E_hat*v     E_hat*v     0 0 0;...
   E_hat*v     E_hat*(1-v) E_hat*v     0 0 0;... 
   E_hat*v     E_hat*v     E_hat*(1-v) 0 0 0;... 
   0           0           0           G 0 0;...
   0           0           0           0 G 0;...
   0           0           0           0 0 G];

s=D*e;

Material_State.e=e;
Material_State.s=s;


end

