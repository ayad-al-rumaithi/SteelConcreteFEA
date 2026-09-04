function [J_det, B]=truss2(V,r)

% Shape functions used:
% N(1)=1/2*(1-r);
% N(2)=1/2*(1+r);

% Use these coefficients if you want to use 1 integration point:  
% r=[0];
% w=[2];

x=V(:,1); 
y=V(:,2); 
z=V(:,3);

le=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2+(z(2)-z(1))^2);
l=(x(2)-x(1))/le;
m=(y(2)-y(1))/le;
n=(z(2)-z(1))/le;

J=le/2;
J_det=J;

B1=[-1/2 1/2]/J;

T=[l m n 0 0 0;...
   0 0 0 l m n];

B=B1*T;
    
end
