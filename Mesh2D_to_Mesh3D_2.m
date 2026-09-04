function [Nodes3D,Prisms,Bricks] = Mesh2D_to_Mesh3D_2(Nodes,Triangles,Quads,zz)

Triangles(:,4)=-10e10;
Mesh2D=[Triangles; Quads];

n=size(Nodes,1);
Nz=length(zz); 

Nodes3D=[];
Bricks=[];

for i=1:1:Nz-1
    
Nodes3D=[Nodes3D; [Nodes zz(i)*ones(n,1)]];   

Bricks=[Bricks; [Mesh2D+(i-1)*n Mesh2D+i*n]];
    
end

Nodes3D=[Nodes3D; [Nodes zz(Nz)*ones(n,1)]];   

A=find(sum(Bricks,2)<0);
Prisms=Bricks(A,[1 2 3 5 6 7]);
Bricks(A,:)=[];

end

