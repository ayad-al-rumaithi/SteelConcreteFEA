function [J_det, B]=spring2(V,Direction,r)

l=Direction(1);
m=Direction(2);
n=Direction(3);

J_det=1;
B=[-l -m -n l m n];
    
end
