function [J_det, B]=wdg6(V,r,s,t)

x=V(:,1); 
y=V(:,2); 
z=V(:,3);

% N(1)=1/2*(1-r-s)*(1-t);
% N(2)=1/2*r*(1-t);
% N(3)=1/2*s*(1-t);
% N(4)=1/2*(1-r-s)*(1+t);
% N(5)=1/2*r*(1+t);
% N(6)=1/2*s*(1+t);

d_N_r=[ t/2 - 1/2; 1/2 - t/2; 0; - t/2 - 1/2; t/2 + 1/2; 0];
d_N_s=[ t/2 - 1/2; 0; 1/2 - t/2; - t/2 - 1/2; 0; t/2 + 1/2];
d_N_t=[ r/2 + s/2 - 1/2; -r/2; -s/2; 1/2 - s/2 - r/2; r/2; s/2];

J=[x'*d_N_r y'*d_N_r z'*d_N_r;...
   x'*d_N_s y'*d_N_s z'*d_N_s;...  
   x'*d_N_t y'*d_N_t z'*d_N_t];
   
J_det=det(J);

d_N_x=zeros(6,1);
d_N_y=zeros(6,1);
d_N_z=zeros(6,1);
for i=1:1:6                                          
NN=J\[d_N_r(i) ; d_N_s(i) ; d_N_t(i)];
d_N_x(i)=NN(1);
d_N_y(i)=NN(2);
d_N_z(i)=NN(3);
end

B=[ d_N_x(1) 0          0  d_N_x(2)    0    0    d_N_x(3)    0    0    d_N_x(4)    0    0    d_N_x(5)    0    0    d_N_x(6)    0    0 ;...
    0       d_N_y(1)    0  0    d_N_y(2)    0    0    d_N_y(3)    0    0    d_N_y(4)    0    0    d_N_y(5)    0    0    d_N_y(6)    0 ;...
    0       0    d_N_z(1)  0    0    d_N_z(2)    0    0    d_N_z(3)    0    0    d_N_z(4)    0    0    d_N_z(5)    0    0    d_N_z(6) ;...
    d_N_y(1) d_N_x(1) 0    d_N_y(2) d_N_x(2) 0   d_N_y(3) d_N_x(3) 0   d_N_y(4) d_N_x(4) 0   d_N_y(5) d_N_x(5) 0   d_N_y(6) d_N_x(6) 0;...
    0 d_N_z(1) d_N_y(1)    0 d_N_z(2) d_N_y(2)   0 d_N_z(3) d_N_y(3)   0 d_N_z(4) d_N_y(4)   0 d_N_z(5) d_N_y(5)   0 d_N_z(6) d_N_y(6);...
    d_N_z(1) 0 d_N_x(1)    d_N_z(2) 0 d_N_x(2)   d_N_z(3) 0 d_N_x(3)   d_N_z(4) 0 d_N_x(4)   d_N_z(5) 0 d_N_x(5)   d_N_z(6) 0 d_N_x(6)];
    
end
