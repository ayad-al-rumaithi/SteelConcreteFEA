function Plot_Results(Model,Show_Deformed,Mag_Factor,Result_type,Materials)

% Result_type:
%   0  - Geometry only
%   1  - Strain component e(1)
%   2  - Strain component e(2)
%   3  - Strain component e(3)
%   4  - Strain component e(4)
%   5  - Strain component e(5)
%   6  - Strain component e(6)
%   7  - Stress component s(1)
%   8  - Stress component s(2)
%   9  - Stress component s(3)
%  10  - Stress component s(4)
%  11  - Stress component s(5)
%  12  - Stress component s(6)
%  13  - Plastic internal variable k
%  14  - Ductile damage history k_D
%  15  - Ductile damage variable d
%  16  - Concrete tensile damage history k_Dt
%  17  - Concrete compressive damage history k_Dc
%  18  - Concrete tensile damage variable d_t
%  19  - Concrete compressive damage variable d_c
%
% Optional Materials argument: vector of material indices to display.
% Omit Materials or pass [] to display all materials.

if nargin<5
    Materials=[];
end

NNode=length(Model.Node);

for i=1:1:NNode

X(i,1)=Model.Node{i}.x;
Y(i,1)=Model.Node{i}.y;
Z(i,1)=Model.Node{i}.z;

Ux(i,1)=Model.Node{i}.Ux;
Uy(i,1)=Model.Node{i}.Uy;
Uz(i,1)=Model.Node{i}.Uz;
    
end

figure;

if Show_Deformed==1
X=X+Ux*Mag_Factor;
Y=Y+Uy*Mag_Factor;
Z=Z+Uz*Mag_Factor;
end

NElement=length(Model.Element);
c=0;
j_Prisms=0;
j_Tetrahedrons=0;
j_Bricks=0;
j_Lines=0;

for i=1:1:NElement
    
if ~isempty(Materials) && ~ismember(Model.Element{i}.Material,Materials)
    continue;
end

c=c+1;
Result(c,1)=0;
e=zeros(6,1); s=zeros(6,1); k=0; k_D=0; d=0; k_Dt=0; k_Dc=0; d_t=0; d_c=0;

No_Int_Points=Model.Element{i}.No_Int_Points;
for j=1:1:No_Int_Points
 
e=e+Model.Element{i}.Integration_Points{j}.e/No_Int_Points;    
s=s+Model.Element{i}.Integration_Points{j}.s/No_Int_Points;    

if strcmp(Model.Material{Model.Element{i}.Material}.Type,'Elastic Plastic Model') || strcmp(Model.Material{Model.Element{i}.Material}.Type,'Elastic Plastic Model 1D')   
k=k+Model.Element{i}.Integration_Points{j}.k/No_Int_Points;                
end

if strcmp(Model.Material{Model.Element{i}.Material}.Type,'Damage Plasticity Model')
k=k+Model.Element{i}.Integration_Points{j}.k/No_Int_Points;
k_Dt=k_Dt+Model.Element{i}.Integration_Points{j}.k_Dt/No_Int_Points;
k_Dc=k_Dc+Model.Element{i}.Integration_Points{j}.k_Dc/No_Int_Points;
d_t=d_t+Model.Element{i}.Integration_Points{j}.d_t/No_Int_Points;
d_c=d_c+Model.Element{i}.Integration_Points{j}.d_c/No_Int_Points;
end

if strcmp(Model.Material{Model.Element{i}.Material}.Type,'Elastic Plastic Ductile Damage Model')    
k=k+Model.Element{i}.Integration_Points{j}.k/No_Int_Points;   
k_D=k_D+Model.Element{i}.Integration_Points{j}.k_D/No_Int_Points;   
d=d+Model.Element{i}.Integration_Points{j}.d/No_Int_Points;   
end

end

if Result_type==1    
Result(c,1)=e(1);
end
if Result_type==2    
Result(c,1)=e(2);
end
if Result_type==3    
Result(c,1)=e(3);
end
if Result_type==4    
Result(c,1)=e(4);
end
if Result_type==5    
Result(c,1)=e(5);
end
if Result_type==6    
Result(c,1)=e(6);
end
if Result_type==7    
Result(c,1)=s(1);
end
if Result_type==8   
Result(c,1)=s(2);
end
if Result_type==9    
Result(c,1)=s(3);
end
if Result_type==10    
Result(c,1)=s(4);
end
if Result_type==11    
Result(c,1)=s(5);
end
if Result_type==12    
Result(c,1)=s(6);
end
if Result_type==13
Result(c,1)=k;
end
if Result_type==14
Result(c,1)=k_D;
end
if Result_type==15
Result(c,1)=d;
end
if Result_type==16
Result(c,1)=k_Dt;
end
if Result_type==17
Result(c,1)=k_Dc;
end
if Result_type==18
Result(c,1)=d_t;
end
if Result_type==19
Result(c,1)=d_c;
end

if strcmp(Model.Element{i,1}.Type,'wdg6')
j_Prisms=j_Prisms+1;    
Prisms(j_Prisms,:)=Model.Element{i}.Nodes;
I_prisms(j_Prisms)=c;
end

if strcmp(Model.Element{i,1}.Type,'tetra4')
j_Tetrahedrons=j_Tetrahedrons+1;   
Tetrahedrons(j_Tetrahedrons,:)=Model.Element{i}.Nodes;
I_Tetrahedrons(j_Tetrahedrons)=c;
end

if strcmp(Model.Element{i,1}.Type,'brick8')...
|| strcmp(Model.Element{i,1}.Type,'brick8 TL')        
j_Bricks=j_Bricks+1;
Bricks(j_Bricks,:)=Model.Element{i}.Nodes;
I_Bricks(j_Bricks)=c;
end   

if strcmp(Model.Element{i,1}.Type,'truss2')
j_Lines=j_Lines+1;
Lines(j_Lines,:)=Model.Element{i}.Nodes;
I_Lines(j_Lines)=c;
end 

end
%--------------------------------------------------------------------------
if j_Prisms>0
Triangles=[];
Rectangles=[];
Triangles=[Triangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,3)]];
Triangles=[Triangles ; [Prisms(:,4) Prisms(:,5) Prisms(:,6)]];
Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,2) Prisms(:,5)  Prisms(:,4)]];
Rectangles=[Rectangles ; [Prisms(:,2) Prisms(:,3) Prisms(:,6)  Prisms(:,5)]];
Rectangles=[Rectangles ; [Prisms(:,1) Prisms(:,3) Prisms(:,6)  Prisms(:,4)]];
if Result_type==0
patch('Faces',Triangles,'Vertices',[X Y Z],'FaceColor','green','EdgeAlpha',0.05);   
patch('Faces',Rectangles,'Vertices',[X Y Z],'FaceColor','green','EdgeAlpha',0.05);
else
ResultT=[Result(I_prisms) ; Result(I_prisms)];    
ResultR=[Result(I_prisms) ; Result(I_prisms) ; Result(I_prisms)];
patch('Faces',Triangles,'Vertices',[X Y Z],'FaceVertexCData',ResultT,'FaceColor','flat','EdgeAlpha',0.05);   
patch('Faces',Rectangles,'Vertices',[X Y Z],'FaceVertexCData',ResultR,'FaceColor','flat','EdgeAlpha',0.05);
colorbar
end    
end    
%--------------------------------------------------------------------------
if j_Tetrahedrons>0
Triangles=[];
Triangles=[Triangles ; [Tetrahedrons(:,1) Tetrahedrons(:,2) Tetrahedrons(:,3)]]; 
Triangles=[Triangles ; [Tetrahedrons(:,1) Tetrahedrons(:,2) Tetrahedrons(:,4)]]; 
Triangles=[Triangles ; [Tetrahedrons(:,2) Tetrahedrons(:,3) Tetrahedrons(:,4)]]; 
Triangles=[Triangles ; [Tetrahedrons(:,3) Tetrahedrons(:,1) Tetrahedrons(:,4)]]; 
if Result_type==0
patch('Faces',Triangles,'Vertices',[X Y Z],'FaceColor','green','EdgeAlpha',0.05);    
else
ResultT=[Result(I_Tetrahedrons) ; Result(I_Tetrahedrons) ; Result(I_Tetrahedrons) ; Result(I_Tetrahedrons)];   
patch('Faces',Triangles,'Vertices',[X Y Z],'FaceVertexCData',ResultT,'FaceColor','flat','EdgeAlpha',0.05);    
colorbar
end    
end
%--------------------------------------------------------------------------
if j_Bricks>0
Rectangles=[];
Rectangles=[Rectangles ; [Bricks(:,1) Bricks(:,2) Bricks(:,3)  Bricks(:,4)]];
Rectangles=[Rectangles ; [Bricks(:,5) Bricks(:,6) Bricks(:,7)  Bricks(:,8)]];
Rectangles=[Rectangles ; [Bricks(:,1) Bricks(:,2) Bricks(:,6)  Bricks(:,5)]];
Rectangles=[Rectangles ; [Bricks(:,2) Bricks(:,3) Bricks(:,7)  Bricks(:,6)]];
Rectangles=[Rectangles ; [Bricks(:,3) Bricks(:,4) Bricks(:,8)  Bricks(:,7)]];
Rectangles=[Rectangles ; [Bricks(:,4) Bricks(:,1) Bricks(:,5)  Bricks(:,8)]];
if Result_type==0
patch('Faces',Rectangles,'Vertices',[X Y Z],'FaceColor','green','EdgeAlpha',0.05);   
else
ResultR=[Result(I_Bricks) ; Result(I_Bricks) ; Result(I_Bricks) ; Result(I_Bricks) ; Result(I_Bricks) ; Result(I_Bricks)];
patch('Faces',Rectangles,'Vertices',[X Y Z],'FaceVertexCData',ResultR,'FaceColor','flat','EdgeAlpha',0.05);    
colorbar
end    
end   
%--------------------------------------------------------------------------
if j_Lines>0
if Result_type==0
for j = 1:j_Lines
line(X(Lines(j,:)),Y(Lines(j,:)),Z(Lines(j,:)),'color','red','LineWidth',1);  
end
else
ResultL= Result(I_Lines);
cmap = colormap;
if max(Result)==min(Result)
    Color_Index=round((size(cmap,1)+1)/2)*ones(size(ResultL));
else
    Color_Index=round(1+(size(cmap,1)-1)* ...
        (ResultL-min(Result))/(max(Result)-min(Result)));
end
for j = 1:j_Lines
line(X(Lines(j,:)),Y(Lines(j,:)),Z(Lines(j,:)),'color',cmap(Color_Index(j),:),'LineWidth',1);   
end  
colorbar  
if max(Result)==min(Result)
    caxis([Result(1)-1,Result(1)+1]);
else
    caxis([min(Result),max(Result)]);
end
end
end
%--------------------------------------------------------------------------
daspect([1 1 1]);
view(52.5,30);

end