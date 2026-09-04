function Model=Model_Set_Initial_Conditions(Model)

disp('-Set Initial Conditions');

Initial_Elastic_Model_State.e=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Model_State.s=[0; 0; 0; 0; 0; 0];

Initial_Elastic_Plastic_Model_State.e=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Plastic_Model_State.s=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Plastic_Model_State.k=0;

Initial_Damage_Plasticity_Model_State.e=[0; 0; 0; 0; 0; 0];
Initial_Damage_Plasticity_Model_State.s=[0; 0; 0; 0; 0; 0];
Initial_Damage_Plasticity_Model_State.s_eff=[0; 0; 0; 0; 0; 0];
Initial_Damage_Plasticity_Model_State.k=0;
Initial_Damage_Plasticity_Model_State.k_D=0;
Initial_Damage_Plasticity_Model_State.d=0;

Initial_Elastic_Plastic_Ductile_Damage_Model.e=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Plastic_Ductile_Damage_Model.s=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Plastic_Ductile_Damage_Model.s_eff=[0; 0; 0; 0; 0; 0];
Initial_Elastic_Plastic_Ductile_Damage_Model.k=0;
Initial_Elastic_Plastic_Ductile_Damage_Model.k_D=0;
Initial_Elastic_Plastic_Ductile_Damage_Model.d=0;

Initial_Elastic_Model_1D_State.e=0;
Initial_Elastic_Model_1D_State.s=0;

Initial_Elastic_Plastic_Model_1D_State.e=0;
Initial_Elastic_Plastic_Model_1D_State.s=0;
Initial_Elastic_Plastic_Model_1D_State.k=0;

Initial_Traction_Seperation_Model_1D_State.e=0;
Initial_Traction_Seperation_Model_1D_State.s=0;
Initial_Traction_Seperation_Model_1D_State.k=0;

for i=1:1:length(Model.Node)
Model.Node{i}.Ux=0; 
Model.Node{i}.Uy=0; 
Model.Node{i}.Uz=0; 
end

for i=1:1:length(Model.Element)

m=Model.Element{i,1}.Material;

if strcmp(Model.Material{m}.Type,'Elastic Model')
Material_State=Initial_Elastic_Model_State;
end
if strcmp(Model.Material{m}.Type,'Elastic Plastic Model')
Material_State=Initial_Elastic_Plastic_Model_State;
end
if strcmp(Model.Material{m}.Type,'Damage Plasticity Model')
Material_State=Initial_Damage_Plasticity_Model_State;
end
if strcmp(Model.Material{m}.Type,'Elastic Plastic Ductile Damage Model')
Material_State=Initial_Elastic_Plastic_Ductile_Damage_Model;
end
if strcmp(Model.Material{m}.Type,'Elastic Model 1D')
Material_State=Initial_Elastic_Model_1D_State;
end
if strcmp(Model.Material{m}.Type,'Elastic Plastic Model 1D')
Material_State=Initial_Elastic_Plastic_Model_1D_State;
end
if strcmp(Model.Material{m}.Type,'Traction Seperation Model 1D')
Material_State=Initial_Traction_Seperation_Model_1D_State;
end

for j=1:1:Model.Element{i,1}.No_Int_Points
Model.Element{i,1}.Integration_Points{j}=Material_State;
end

end

end