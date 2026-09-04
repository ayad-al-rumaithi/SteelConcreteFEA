function [U_f]=Solve_Linear_System(K_ff,F_eff)

disp('-Solve Linear System');

U_f=full(K_ff\sparse(F_eff));

end