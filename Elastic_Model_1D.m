function [Material_State,D]=Elastic_Model_1D(Material,Material_State,e)

E=Material.E;

D=E;
s=D*e;

Material_State.e=e;
Material_State.s=s;

end
