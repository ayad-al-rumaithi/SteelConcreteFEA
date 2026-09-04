function [Material_State,D]=Damage_Plasticity_Model(Material,Material_State,e)

% STAGE 1:
% Unified Abaqus Concrete Damaged Plasticity (Lubliner-Lee-Fenves)
% yield surface with Abaqus-type hyperbolic non-associated flow potential.
%
% Stage-1 scope:
%   1) Rankine + Drucker-Prager multi-surface formulation is replaced by
%      ONE yield surface.
%   2) One cumulative plastic multiplier/internal variable k is used.
%   3) Existing hardening/softening and damage equations are deliberately
%      NOT replaced yet.
%   4) Existing function input/output interface is retained.
%      For backward compatibility, Material_State.k_DP and k_RK are both
%      kept and synchronized to the single k.
%   5) Additional CDP parameters are constants below and are NOT added to
%      the Material input structure.
%
% Important:
%   The present stage keeps the original scalar damage variable d and its
%   evolution exactly in the final part of the original routine. The
%   MC2010 compression/tension laws and two damage variables are Stage 2+.

% Additional unified-CDP parameters (assumed constants for Stage 1)
% -------------------------------------------------------------------------
psi_deg = 35.0;       % Dilation angle, degrees
ecc      = 0.10;      % Abaqus meridional flow-potential eccentricity
Kc       = 2.0/3.0;   % Deviatoric-plane parameter, 0.5 < Kc <= 1.0

% Numerical safeguards
tol_q    = 1.0e-12;
tol_F   = 1.0e-10;
tol_R   = 1.0e-10;
maxIter  = 100;

% Material properties
% --------------------
E   = Material.E;
v   = Material.v;
f_t = Material.f_t;
g_f = Material.g_f;
f_c = Material.f_c;
f_c2= Material.f_c2;

% The original code used these quantities for its DP hardening surface.
% They are retained because the original material interface is retained.
beta_old = sqrt(3)*(f_c2-f_c)/(2*f_c2-f_c);
Hp       = f_c2*f_c/(sqrt(3)*(2*f_c2-f_c)); %#ok<NASGU>

% Previous increment
% -------------------
e_i = Material_State.e;
s_i = Material_State.s_eff;       % effective stress used by return mapping
k_D = Material_State.k_D;

% Backward-compatible storage:
% one unified internal variable k is represented by both old fields.
if isfield(Material_State,'k_DP')
    k_i = Material_State.k_DP;
elseif isfield(Material_State,'k_RK')
    k_i = Material_State.k_RK;
else
    k_i = 0;
end

% Elastic constitutive matrix
% ----------------------------
E_hat = E/((1-2*v)*(1+v));
G     = 0.5*E/(1+v);

D_e = [E_hat*(1-v) E_hat*v     E_hat*v     0 0 0;...
       E_hat*v     E_hat*(1-v) E_hat*v     0 0 0;...
       E_hat*v     E_hat*v     E_hat*(1-v) 0 0 0;...
       0           0           0           G 0 0;...
       0           0           0           0 G 0;...
       0           0           0           0 0 G];

% Elastic predictor
% ------------------
d_e = e-e_i;
s_e = s_i + D_e*d_e;

% Abaqus CDP parameters
% ----------------------
rb0_rc0  = f_c2/f_c;    % sigma_b0 / sigma_c0
alpha = (rb0_rc0-1)/(2*rb0_rc0-1);
gamma = 3*(1-Kc)/(2*Kc-1);
psi   = psi_deg*pi/180;

% Stage 1: strengths are kept constant. Their evolution will be replaced
% later by the MC2010 compression/tension laws.
sigma_t = f_t;
sigma_c = f_c;

% Abaqus beta parameter
beta = (sigma_c/sigma_t)*(1-alpha) - (1+alpha);

% Trial yield check
% ------------------
[sigma_max_e,~,~,~,~] = CDP_Invariants(s_e,tol_q);

[F_e,~,~,~] = CDP_Yield(s_e,sigma_t,sigma_c,alpha,beta,gamma,tol_q);

if F_e <= tol_F
    % Elastic effective stress
    s_eff = s_e;
    D_ep  = D_e;
    d     = 1-exp(-f_t*k_D/g_f);
    s     = (1-d)*s_eff;
    D     = (1-d)*D_ep;

else
    % Return mapping: ONE yield surface + ONE plastic multiplier
    % ---------------------------------------------------------------
    s = s_e;
    k = k_i;

    converged = false;

    for it = 1:maxIter

        [F, dF_ds, d2F_ds2, ~] = CDP_Yield(s,sigma_t,sigma_c,...
                                             alpha,beta,gamma,tol_q);

        [~, dG_ds, d2G_ds2, ~] = CDP_Flow(s,sigma_t,psi,ecc,tol_q);

        % Residuals:
        %   r_sigma = De^{-1}(s-se) + (k-ki)*G_,s = 0
        %   r_F     = F(s) = 0
        dk = k-k_i;

        r1 = D_e\(s-s_e) + dk*dG_ds;
        r2 = F;
        r  = [r1;r2];

        if norm(r1) <= tol_R*max(1,f_t) && ...
           abs(r2) <= tol_F*max(1,f_t)
            converged = true;
            break
        end

        % Exact-form Newton Jacobian for the chosen constitutive equations:
        %
        % dr_sigma/ds = De^{-1} + dk*d2G/ds2
        % dr_sigma/dk = dG/ds
        % dr_F/ds     = dF/ds
        % dr_F/dk     = 0  (Stage 1: no hardening evolution yet)
        J11 = inv(D_e) + dk*d2G_ds2;
        J12 = dG_ds;
        J21 = dF_ds';
        J22 = 0;
        J   = [J11 J12; J21 J22];

        delta = -J\r;

        ds = delta(1:6);
        dk_new = delta(7);

        % Prevent a Newton step from making the cumulative multiplier
        % negative. Backtracking is used instead of the old multi-surface
        % restart logic.
        if k + dk_new < k_i
            scale = 1.0;
            while k + scale*dk_new < k_i
                scale = 0.5*scale;
                if scale < 1e-8
                    break
                end
            end
            ds     = scale*ds;
            dk_new = scale*dk_new;
        end

        % Simple residual-based backtracking for robustness.
        base_norm = norm(r);
        scale = 1.0;
        accepted = false;

        while scale >= 1e-6
            s_try = s + scale*ds;
            k_try = k + scale*dk_new;

            [F_try,~,~,~] = CDP_Yield(s_try,sigma_t,sigma_c,...
                                       alpha,beta,gamma,tol_q);
            dk_try = k_try-k_i;
            [~,dG_try,~,~] = CDP_Flow(s_try,sigma_t,psi,ecc,tol_q);
            r1_try = D_e\(s_try-s_e) + dk_try*dG_try;
            r_try  = [r1_try;F_try];

            if norm(r_try) < base_norm
                s = s_try;
                k = k_try;
                accepted = true;
                break
            end
            scale = 0.5*scale;
        end

        if ~accepted
            % Accept the Newton step as a last resort; convergence check
            % below will catch a failed iteration.
            s = s + ds;
            k = max(k_i,k+dk_new);
        end
    end

    if ~converged
        error('Damage_Plasticity_Model:ReturnMapping',...
              'Unified CDP return mapping did not converge in %d iterations.',maxIter);
    end

    % Consistent elastoplastic tangent for Stage 1
    % ------------------------------------------------
    dk = k-k_i;

    [F,dF_ds,d2F_ds2,~] = CDP_Yield(s,sigma_t,sigma_c,...
                                     alpha,beta,gamma,tol_q); %#ok<ASGLU>
    [~,dG_ds,d2G_ds2,~] = CDP_Flow(s,sigma_t,psi,ecc,tol_q);

    % H is the stress-space tangent associated with the current return
    % mapping equation.
    Hinv = inv(D_e) + dk*d2G_ds2;
    H    = inv(Hinv);

    % Yield consistency uses F_,s while plastic flow uses G_,s.
    denom = dF_ds'*H*dG_ds;

    if abs(denom) < 1e-14
        error('Damage_Plasticity_Model:SingularTangent',...
              'The unified CDP consistent tangent denominator is too small.');
    end

    D_ep = H - H*dG_ds*(1/denom)*dF_ds'*H;

    % Existing plastic-strain and damage part
    % -----------------------------------------
    e_p = e-D_e\s;
    d_e_p_d_e = inv(D_e)*(D_e-D_ep);

    % Principal Plastic Strain
    [Tp,Ep] = eig([e_p(1) e_p(4)/2 e_p(6)/2;...
                   e_p(4)/2 e_p(2) e_p(5)/2;...
                   e_p(6)/2 e_p(5)/2 e_p(3)]);
    Ep=diag(Ep);
    [Ep,I]=sort(Ep,'descend');
    Tp=Tp(:,I);

    TTp=[Tp(1,1)^2 Tp(2,1)^2 Tp(3,1)^2 Tp(1,1)*Tp(2,1) Tp(2,1)*Tp(3,1) Tp(3,1)*Tp(1,1);...
         Tp(1,2)^2 Tp(2,2)^2 Tp(3,2)^2 Tp(1,2)*Tp(2,2) Tp(2,2)*Tp(3,2) Tp(3,2)*Tp(1,2);...
         Tp(1,3)^2 Tp(2,3)^2 Tp(3,3)^2 Tp(1,3)*Tp(2,3) Tp(2,3)*Tp(3,3) Tp(3,3)*Tp(1,3)];

    Ep_eqv=sqrt(Ep(1)^2+Ep(2)^2+Ep(3)^2);

    if Ep_eqv > 0
        d_Ep_eqv_d_Ep=[Ep(1)/(Ep_eqv);...
                       Ep(2)/(Ep_eqv);...
                       Ep(3)/(Ep_eqv)];
    else
        d_Ep_eqv_d_Ep=zeros(3,1);
    end

    d_Ep_eqv_d_ep=TTp'*d_Ep_eqv_d_Ep;
    d_Ep_eqv_d_e=d_e_p_d_e'*d_Ep_eqv_d_ep;

    % Loading and Unloading Conditions
    if Ep_eqv>k_D
        k_D=Ep_eqv;
        d_k_D_d_Ep_eqv=1;
    else
        d_k_D_d_Ep_eqv=0;
    end

    d_k_D_d_e=d_k_D_d_Ep_eqv*d_Ep_eqv_d_e;

    % Existing damage law retained unchanged
    d=1-exp(-f_t*k_D/g_f);
    d_d_d_k_D=f_t/g_f*exp(-f_t*k_D/g_f);
    d_d_d_e=d_d_d_k_D*d_k_D_d_e;

    % Stress and tangent
    s_eff=s;
    s=(1-d)*s_eff;
    D=(1-d)*D_ep-s_eff*d_d_d_e';
end

% Final material state
% ---------------------
Material_State.e     = e;
Material_State.s     = s;
Material_State.s_eff = s_eff;

% Keep the old output fields, but synchronize both to the ONE Stage-1 k.
Material_State.k_RK = k_i;
Material_State.k_DP = k_i;

if exist('k','var')
    Material_State.k_RK = k;
    Material_State.k_DP = k;
end

Material_State.k_D = k_D;
Material_State.d   = d;

end


% ========================================================================
function [F,dF_ds,d2F_ds2,aux] = CDP_Yield(s,sigma_t,sigma_c,alpha,beta,gamma,tol_q)
% Abaqus CDP Lubliner-Lee-Fenves yield surface.
%
% F = 1/(1-alpha) * [ q - 3*alpha*p
%                     + beta*<sigma_max>
%                     - gamma*<-sigma_max> ] - sigma_c
%
% The derivatives are evaluated directly in Cartesian/Voigt stress space.
% d2F/ds2 is obtained by a symmetric central difference of the analytical
% first derivative. This retains the correct dependence of the principal
% stress on the stress tensor and avoids the incomplete T'*H*T operation
% used in the original routine.

[p,q,sdev,sigma_max,Pmax] = CDP_Invariants(s,tol_q);

pos = max(sigma_max,0);
neg = max(-sigma_max,0);

F = (q - 3*alpha*p + beta*pos - gamma*neg)/(1-alpha) - sigma_c;

% Derivative of p
dp = [-1/3;-1/3;-1/3;0;0;0];

% Derivative of q
if q > tol_q
    dq = [3/(2*q)*sdev(1);...
          3/(2*q)*sdev(2);...
          3/(2*q)*sdev(3);...
          3/q*sdev(4);...
          3/q*sdev(5);...
          3/q*sdev(6)];
else
    dq = zeros(6,1);
end

% d sigma_max / d stress = principal projector of maximum eigenvalue
dsm = [Pmax(1,1);Pmax(2,2);Pmax(3,3);...
       2*Pmax(1,2);2*Pmax(2,3);2*Pmax(3,1)];

% At sigma_max=0 the Macaulay bracket is nondifferentiable.
% A centered numerical derivative is preferable to an arbitrary branch.
if abs(sigma_max) <= 10*tol_q
    dbr_pos = 0.5*dsm;
    dbr_neg = -0.5*dsm;
elseif sigma_max > 0
    dbr_pos = dsm;
    dbr_neg = zeros(6,1);
else
    dbr_pos = zeros(6,1);
    dbr_neg = -dsm;
end

dF_ds = (dq - 3*alpha*dp + beta*dbr_pos - gamma*dbr_neg)/(1-alpha);

% Numerically differentiate the COMPLETE first derivative in Cartesian
% stress space. This includes the dependence of principal directions.
h = max(1e-7*max([sigma_t; sigma_c; q; 1]),1e-10);
d2F_ds2 = zeros(6,6);

for a=1:6
    ds = zeros(6,1);
    ds(a)=h;
    [~,gp,~,~] = CDP_Yield_FirstDerivative(s+ds,sigma_t,sigma_c,...
                                            alpha,beta,gamma,tol_q);
    [~,gm,~,~] = CDP_Yield_FirstDerivative(s-ds,sigma_t,sigma_c,...
                                            alpha,beta,gamma,tol_q);
    d2F_ds2(:,a)=(gp-gm)/(2*h);
end

% Enforce the expected symmetry of the Hessian up to numerical error.
d2F_ds2=0.5*(d2F_ds2+d2F_ds2');

aux.p=p; aux.q=q; aux.sdev=sdev; aux.sigma_max=sigma_max;
end


function [F,dF_ds,d2F_ds2,aux] = CDP_Yield_FirstDerivative(s,sigma_t,sigma_c,alpha,beta,gamma,tol_q)
[p,q,sdev,sigma_max, Pmax] = CDP_Invariants(s,tol_q);

F=(q-3*alpha*p+beta*max(sigma_max,0)-...
   gamma*max(-sigma_max,0))/(1-alpha)-sigma_c;

dp=[-1/3;-1/3;-1/3;0;0;0];

if q>tol_q
    dq=[3/(2*q)*sdev(1);...
        3/(2*q)*sdev(2);...
        3/(2*q)*sdev(3);...
        3/q*sdev(4);...
        3/q*sdev(5);...
        3/q*sdev(6)];
else
    dq=zeros(6,1);
end

dsm=[Pmax(1,1);Pmax(2,2);Pmax(3,3);...
     2*Pmax(1,2);2*Pmax(2,3);2*Pmax(3,1)];

if abs(sigma_max)<=10*tol_q
    dbr_pos=0.5*dsm;
    dbr_neg=-0.5*dsm;
elseif sigma_max>0
    dbr_pos=dsm; dbr_neg=zeros(6,1);
else
    dbr_pos=zeros(6,1); dbr_neg=-dsm;
end

dF_ds=(dq-3*alpha*dp+beta*dbr_pos-gamma*dbr_neg)/(1-alpha);
d2F_ds2=zeros(6,6);
aux=[];
end


function [G,dG_ds,d2G_ds2,aux] = CDP_Flow(s,sigma_t,psi,ecc,tol_q)
% Abaqus CDP hyperbolic flow potential:
% G = sqrt((ecc*sigma_t*tan(psi))^2 + q^2) - p*tan(psi)

[p,q,sdev,~,~]=CDP_Invariants(s,tol_q);

A=ecc*sigma_t*tan(psi);
R=sqrt(A^2+q^2);

G=R-p*tan(psi);

dp=[-1/3;-1/3;-1/3;0;0;0];

if q>tol_q
    dq=[3/(2*q)*sdev(1);...
        3/(2*q)*sdev(2);...
        3/(2*q)*sdev(3);...
        3/q*sdev(4);...
        3/q*sdev(5);...
        3/q*sdev(6)];
else
    dq=zeros(6,1);
end

dG_ds=(q/R)*dq-tan(psi)*dp;

% Exact Hessian of G through the Hessian of q:
% d2G = (A^2/R^3) dq*dq' + (q/R)*d2q
%
% q=sqrt(3/2 s:s), hence in tensor/Voigt representation:
% d2q is the deviatoric metric term projected consistently into the
% engineering-shear Voigt representation.
if q>tol_q
    % Numerical Hessian of q gives a robust implementation for the
    % engineering-shear convention used by the original MATLAB code.
    h=max(1e-7*max([sigma_t;q;1]),1e-10);
    d2q=zeros(6,6);
    for a=1:6
        ds=zeros(6,1); ds(a)=h;
        [~,gp]=Q_FirstDerivative(s+ds,tol_q);
        [~,gm]=Q_FirstDerivative(s-ds,tol_q);
        d2q(:,a)=(gp-gm)/(2*h);
    end
    d2q=0.5*(d2q+d2q');
    d2G_ds2=(A^2/R^3)*(dq*dq')+(q/R)*d2q;
else
    d2G_ds2=zeros(6,6);
end

aux.p=p; aux.q=q; aux.A=A;
end


function [q,dq]=Q_FirstDerivative(s,tol_q)
[p,q,sdev,~,~]=CDP_Invariants(s,tol_q); %#ok<ASGLU>
if q>tol_q
    dq=[3/(2*q)*sdev(1);...
        3/(2*q)*sdev(2);...
        3/(2*q)*sdev(3);...
        3/q*sdev(4);...
        3/q*sdev(5);...
        3/q*sdev(6)];
else
    dq=zeros(6,1);
end
end


function [p,q,sdev,sigma_max,Pmax]=CDP_Invariants(s,tol_q)
% Voigt order: [11 22 33 12 23 31]
Smat=[s(1) s(4) s(6);...
      s(4) s(2) s(5);...
      s(6) s(5) s(3)];

[V,L]=eig(Smat);
lam=diag(L);
[lam,idx]=sort(lam,'descend');
V=V(:,idx);

sigma_max=lam(1);
Pmax=V(:,1)*V(:,1)';

I1=trace(Smat);
p=-I1/3;

sdev=Smat+ p*eye(3);

J2tensor=0.5*sum(sum(sdev.*sdev));
q=sqrt(max(0,3*J2tensor));

if q<tol_q
    sdev=[sdev(1,1);sdev(2,2);sdev(3,3);...
          sdev(1,2);sdev(2,3);sdev(3,1)];
else
    sdev=[sdev(1,1);sdev(2,2);sdev(3,3);...
          sdev(1,2);sdev(2,3);sdev(3,1)];
end
end
