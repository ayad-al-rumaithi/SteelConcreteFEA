function Force=Get_Force(F_ex)

c=0;
Force=[];
for i=1:1:length(F_ex);
if F_ex(i)~=0
c=c+1;    
Force{c}.DOF=i; Force{c}.Magnitude=F_ex(i);
end
end

end