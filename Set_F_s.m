function Model=Set_F_s(Model,F_s)
disp('-Update Model Constraint Forces');
for i=1:1:length(Model.Constraint)                        
Model.Constraint{i,1}.F=F_s(i);
end

end

