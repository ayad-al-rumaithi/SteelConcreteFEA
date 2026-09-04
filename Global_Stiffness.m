function K=Global_Stiffness(Model)  
disp('-Global Assembly');
NDOF=Model.NDOF;
NElement=length(Model.Element);
mk=0;

for i=1:1:NElement
    
NNodeElement=length(Model.Element{i,1}.Nodes);    

NN=[];
for j=1:1:NNodeElement   
Node_num=Model.Element{i,1}.Nodes(j);                             
NN(3*j-2)=Model.Node{Node_num}.DOFx;  NN(3*j-1)=Model.Node{Node_num}.DOFy;  NN(3*j)=Model.Node{Node_num}.DOFz;    
end  

NDOFElement=length(NN);
M=meshgrid(1:NDOFElement)';M2=M';
ii=M(:);
jj=M2(:);
kk=1:NDOFElement^2;
k=Model.Element{i,1}.K;   
kkk=kk+mk; 
mk=mk+NDOFElement^2;

Ig(kkk)=NN(ii);
Jg(kkk)=NN(jj);
Kg(kkk)=k(:);

end
 
K=sparse(Ig,Jg,Kg,NDOF,NDOF);

end

