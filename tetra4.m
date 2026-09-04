function [J_det, B]=tetra4(V,r,s,t)

% Shape functions used:
% N(1)=1-r-s-t;
% N(2)=r;
% N(3)=s;
% N(4)=t;

% Use these coefficients if you want to use 1 integration point:  
%r=[1/4];
%s=[1/4];
%t=[1/4];
%w=[1/6];

x=V(:,1); 
y=V(:,2); 
z=V(:,3);

d_N_r=[-1; 1; 0; 0];
d_N_s=[-1; 0; 1; 0];
d_N_t=[-1; 0; 0; 1];

J=[x'*d_N_r y'*d_N_r z'*d_N_r;...
   x'*d_N_s y'*d_N_s z'*d_N_s;...  
   x'*d_N_t y'*d_N_t z'*d_N_t];
   
J_det=det(J);

d_N_x=zeros(4,1);
d_N_y=zeros(4,1);
d_N_z=zeros(4,1);
for i=1:1:4                                          
NN=J\[d_N_r(i) ; d_N_s(i) ; d_N_t(i)];
d_N_x(i)=NN(1);
d_N_y(i)=NN(2);
d_N_z(i)=NN(3);
end

B=[ d_N_x(1) 0          0  d_N_x(2)    0    0    d_N_x(3)    0    0    d_N_x(4)    0    0  ;...
    0       d_N_y(1)    0  0    d_N_y(2)    0    0    d_N_y(3)    0    0    d_N_y(4)    0  ;...
    0       0    d_N_z(1)  0    0    d_N_z(2)    0    0    d_N_z(3)    0    0    d_N_z(4)  ;...
    d_N_y(1) d_N_x(1) 0    d_N_y(2) d_N_x(2) 0   d_N_y(3) d_N_x(3) 0   d_N_y(4) d_N_x(4) 0 ;...
    0 d_N_z(1) d_N_y(1)    0 d_N_z(2) d_N_y(2)   0 d_N_z(3) d_N_y(3)   0 d_N_z(4) d_N_y(4) ;...
    d_N_z(1) 0 d_N_x(1)    d_N_z(2) 0 d_N_x(2)   d_N_z(3) 0 d_N_x(3)   d_N_z(4) 0 d_N_x(4)];
    
end

