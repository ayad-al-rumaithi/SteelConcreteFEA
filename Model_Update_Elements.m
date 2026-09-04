function Model=Model_Update_Elements(Model,Modeli) 


for i=1:1:length(Model.Element)
    
V=[]; U=[];

for j=1:1:length(Model.Element{i,1}.Nodes)
    
Node_num=Model.Element{i,1}.Nodes(j);

V(j,1)=Model.Node{Node_num}.x;  V(j,2)=Model.Node{Node_num}.y;  V(j,3)=Model.Node{Node_num}.z;

U(j,1)=Model.Node{Node_num}.Ux; U(j,2)=Model.Node{Node_num}.Uy; U(j,3)=Model.Node{Node_num}.Uz;   
    
end

Element=Modeli.Element{i,1};
Material=Modeli.Material{Element.Material};

[Element]=Element_Analysis(Element,Material,V,U);

Model.Element{i,1}=Element;
  
end

end




