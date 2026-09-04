function [Material_State,D]=Material_Model(Material,Material_State,e)

%Elastic Model
%-------------------------------------

if strcmp(Material.Type,'Elastic Model')   

[Material_State,D]=Elastic_Model(Material,Material_State,e);

end

%Elastic Plastic Model
%-------------------------------------

if strcmp(Material.Type,'Elastic Plastic Model')   

[Material_State,D]=Elastic_Plastic_Model(Material,Material_State,e);
    
end    

%Damage Plasticity Model
%-------------------------------------

if strcmp(Material.Type,'Damage Plasticity Model')   

[Material_State,D]=Damage_Plasticity_Model(Material,Material_State,e);

end

%Elastic Plastic Ductile Damage Model
%-------------------------------------

if strcmp(Material.Type,'Elastic Plastic Ductile Damage Model')   

[Material_State,D]=Elastic_Plastic_Ductile_Damage_Model(Material,Material_State,e);

end

%Elastic Model 1D
%-------------------------------------

if strcmp(Material.Type,'Elastic Model 1D')   

[Material_State,D]=Elastic_Model_1D(Material,Material_State,e);

end

%Elastic Plastic Model 1D
%-------------------------------------

if strcmp(Material.Type,'Elastic Plastic Model 1D')   

[Material_State,D]=Elastic_Plastic_Model_1D(Material,Material_State,e);

end

%Traction Seperation Model 1D
%-------------------------------------

if strcmp(Material.Type,'Traction Seperation Model 1D')   

[Material_State,D]=Traction_Seperation_Model_1D(Material,Material_State,e);

end

end
