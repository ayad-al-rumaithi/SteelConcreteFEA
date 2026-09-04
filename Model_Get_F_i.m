function [F_i]=Model_Get_F_i(Model)

NDOF=Model.NDOF;
F_i=zeros(NDOF,1);

for i=1:1:length(Model.Element)

f=Model.Element{i,1}.f;

if isnan(f); error('Internal Force Value is NaN'); end

NN=[];
for j=1:1:length(Model.Element{i,1}.Nodes)
Node_num=Model.Element{i,1}.Nodes(j);
NN(3*j-2)=Model.Node{Node_num}.DOFx;  NN(3*j-1)=Model.Node{Node_num}.DOFy;  NN(3*j)=Model.Node{Node_num}.DOFz;    
end  

F_i(NN)=F_i(NN)+f;
end

end