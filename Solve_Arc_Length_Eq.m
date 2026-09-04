function [dlambda]=Solve_Arc_Length_Eq(a1,a2,a3,a4,a5,det_K)

disp('-Solve Arc Length Equation');
if (a1==0 && a2==0)
    error('No roots found')
elseif (a1==0 && a2~=0)
    dlambda=-a3/a2;  %Linear Solution
else
    fac=a2*a2-4*a1*a3;
    if fac<0
    dlambda=-a2/(2*a1);  %Global Optimium
    disp('--***Imaginary Roots Found, Global Optimium was Used');
    else
    dlambda1=(-a2+sqrt(fac))/(2*a1);
    dlambda2=(-a2-sqrt(fac))/(2*a1);
    
    DOT1=a4+dlambda1*a5;
    DOT2=a4+dlambda2*a5;
    
    dlambda=dlambda1;
    if DOT2>DOT1; dlambda=dlambda2; end;
    
    if a4==0 && a5==0; 
        if det_K*dlambda1<0 && det_K*dlambda2>0
            dlambda=dlambda2;
        end
    end
        
    end

end


end
