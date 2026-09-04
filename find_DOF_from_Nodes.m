function [DOFx, DOFy, DOFz, Node]=find_DOF_from_Nodes(Nodes3D,x,y,z)

DOFx=[];
DOFy=[];

NNode=size(Nodes3D,1);
c=0;

for i=1:1:NNode
    
if not(isempty(x)); if abs(Nodes3D(i,1)-x)>10^(-6); continue; end; end
if not(isempty(y)); if abs(Nodes3D(i,2)-y)>10^(-6); continue; end; end
if not(isempty(z)); if abs(Nodes3D(i,3)-z)>10^(-6); continue; end; end

c=c+1;

DOFx(c)=3*i-2;
DOFy(c)=3*i-1;
DOFz(c)=3*i;
Node(c)=i;

end

end





