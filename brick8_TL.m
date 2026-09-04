function [J_det, B, G]=brick8_TL(V,UU,r,s,t)

% Shape functions used:
% N(1)=1/8*(1-r)*(1-s)*(1-t);
% N(2)=1/8*(1+r)*(1-s)*(1-t);
% N(3)=1/8*(1+r)*(1+s)*(1-t);
% N(4)=1/8*(1-r)*(1+s)*(1-t);
% N(5)=1/8*(1-r)*(1-s)*(1+t);
% N(6)=1/8*(1+r)*(1-s)*(1+t);
% N(7)=1/8*(1+r)*(1+s)*(1+t);
% N(8)=1/8*(1-r)*(1+s)*(1+t);

% Use these coefficients if you want to use 8 integration points:  
%r=1/sqrt(3)*[-1; 1; 1; -1; -1; 1; 1; -1];
%s=1/sqrt(3)*[-1; -1; 1; 1; -1; -1; 1; 1];
%t=1/sqrt(3)*[-1; -1; -1; -1; 1; 1; 1; 1];
%w=[1; 1; 1; 1; 1; 1; 1; 1];

x=V(:,1); 
y=V(:,2); 
z=V(:,3);

d_N_r=[ -((s - 1)*(t - 1))/8; ((s - 1)*(t - 1))/8; -((s + 1)*(t - 1))/8; ((s + 1)*(t - 1))/8; ((s - 1)*(t + 1))/8; -((s - 1)*(t + 1))/8; ((s + 1)*(t + 1))/8; -((s + 1)*(t + 1))/8];
d_N_s=[ -(r/8 - 1/8)*(t - 1); (r/8 + 1/8)*(t - 1); -(r/8 + 1/8)*(t - 1); (r/8 - 1/8)*(t - 1); (r/8 - 1/8)*(t + 1); -(r/8 + 1/8)*(t + 1); (r/8 + 1/8)*(t + 1); -(r/8 - 1/8)*(t + 1)];
d_N_t=[ -(r/8 - 1/8)*(s - 1); (r/8 + 1/8)*(s - 1); -(r/8 + 1/8)*(s + 1); (r/8 - 1/8)*(s + 1); (r/8 - 1/8)*(s - 1); -(r/8 + 1/8)*(s - 1); (r/8 + 1/8)*(s + 1); -(r/8 - 1/8)*(s + 1)];

J=[x'*d_N_r y'*d_N_r z'*d_N_r;...
   x'*d_N_s y'*d_N_s z'*d_N_s;...  
   x'*d_N_t y'*d_N_t z'*d_N_t];
   
J_det=det(J);

d_N_x=zeros(8,1);
d_N_y=zeros(8,1);
d_N_z=zeros(8,1);
for i=1:1:8                                          
NN=J\[d_N_r(i) ; d_N_s(i) ; d_N_t(i)];
d_N_x(i)=NN(1);
d_N_y(i)=NN(2);
d_N_z(i)=NN(3);
end

H=[1 0 0 0 0 0 0 0 0;...
   0 0 0 0 1 0 0 0 0;...
   0 0 0 0 0 0 0 0 1;...
   0 1 0 1 0 0 0 0 0;...
   0 0 0 0 0 1 0 1 0;...
   0 0 1 0 0 0 1 0 0];...
   
G=[d_N_x(1) 0       0     d_N_x(2) 0       0    d_N_x(3) 0       0    d_N_x(4) 0       0    d_N_x(5) 0       0    d_N_x(6) 0       0    d_N_x(7) 0       0    d_N_x(8) 0       0;...
   d_N_y(1) 0       0     d_N_y(2) 0       0    d_N_y(3) 0       0    d_N_y(4) 0       0    d_N_y(5) 0       0    d_N_y(6) 0       0    d_N_y(7) 0       0    d_N_y(8) 0       0;...
   d_N_z(1) 0       0     d_N_z(2) 0       0    d_N_z(3) 0       0    d_N_z(4) 0       0    d_N_z(5) 0       0    d_N_z(6) 0       0    d_N_z(7) 0       0    d_N_z(8) 0       0;...
   0    d_N_x(1)    0     0    d_N_x(2)    0    0    d_N_x(3)    0    0    d_N_x(4)    0    0    d_N_x(5)    0    0    d_N_x(6)    0    0    d_N_x(7)    0    0    d_N_x(8)    0;...
   0    d_N_y(1)    0     0    d_N_y(2)    0    0    d_N_y(3)    0    0    d_N_y(4)    0    0    d_N_y(5)    0    0    d_N_y(6)    0    0    d_N_y(7)    0    0    d_N_y(8)    0;...
   0    d_N_z(1)    0     0    d_N_z(2)    0    0    d_N_z(3)    0    0    d_N_z(4)    0    0    d_N_z(5)    0    0    d_N_z(6)    0    0    d_N_z(7)    0    0    d_N_z(8)    0;...
   0    0    d_N_x(1)     0    0    d_N_x(2)    0    0    d_N_x(3)    0    0    d_N_x(4)    0    0    d_N_x(5)    0    0    d_N_x(6)    0    0    d_N_x(7)    0    0    d_N_x(8);...
   0    0    d_N_y(1)     0    0    d_N_y(2)    0    0    d_N_y(3)    0    0    d_N_y(4)    0    0    d_N_y(5)    0    0    d_N_y(6)    0    0    d_N_y(7)    0    0    d_N_y(8);...
   0    0    d_N_z(1)     0    0    d_N_z(2)    0    0    d_N_z(3)    0    0    d_N_z(4)    0    0    d_N_z(5)    0    0    d_N_z(6)    0    0    d_N_z(7)    0    0    d_N_z(8)];...

B_o=H*G;

theta=G*UU;
d_u_x=theta(1); d_u_y=theta(2); d_u_z=theta(3); 
d_v_x=theta(4); d_v_y=theta(5); d_v_z=theta(6); 
d_w_x=theta(7); d_w_y=theta(8); d_w_z=theta(9); 

A=[d_u_x  0      0      d_v_x  0      0      d_w_x  0      0    ;...
   0      d_u_y  0      0      d_v_y  0      0      d_w_y  0    ;...
   0      0      d_u_z  0      0      d_v_z  0      0      d_w_z;...
   d_u_y  d_u_x  0      d_v_y  d_v_x  0      d_w_y  d_w_x  0    ;...
   0      d_u_z  d_u_y  0      d_v_z  d_v_y  0      d_w_z  d_w_y;...
   d_u_z  0      d_u_x  d_v_z  0      d_v_x  d_w_z  0      d_w_x];...
    
B_nl=A*G;   
    
B=B_o+B_nl;

end
