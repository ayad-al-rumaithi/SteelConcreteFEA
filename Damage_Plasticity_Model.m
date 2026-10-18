function [Material_State,D]=Damage_Plasticity_Model(Material,Material_State,e)

% STAGE II-A: optimized separate tension/compression damage (intermediate model).
% Based on verified Stage I. Effective-stress CDP return mapping is unchanged.
% Histories use stress-weighted principal PLASTIC INCREMENTS. d_t is the retained
% provisional exponential tensile law, while d_c = 0 until Stage II-C.
% Nominal stress uses spectral tension/compression split of effective stress.
% This produces unilateral compressive stiffness recovery after tensile damage.
% The stress-split law is nonsmooth at zero principal stress; tangent there
% is a centered generalized derivative. No MC2010 laws are implemented yet.
%
% Required initial state fields: e,s,s_eff,k,k_Dt,k_Dc,d_t,d_c.
% Initial values: zeros(6,1) for e,s,s_eff; zero for other fields.
% No backward compatibility with Stage-I k_D/d is provided.
% Inputs/outputs and material fields remain unchanged.
%
% Additional unified-CDP parameters (assumed constants for Stage 1)
% -------------------------------------------------------------------------
psi_deg = 35.0;       % Dilation angle, degrees
ecc      = 0.10;      % Abaqus meridional flow-potential eccentricity
Kc       = 2.0/3.0;   % Deviatoric-plane parameter, 0.5 < Kc <= 1.0
rb0_rc0 = 1.16;       % Biaxial/uniaxial compressive strength ratio

% Numerical safeguards
% ---------------------
tol_q    = 1.0e-12;
tol_F   = 1.0e-12;
tol_R   = 1.0e-12;
maxIter  = 100;

% Material properties
% --------------------
E   = Material.E;
v   = Material.v;
f_t = Material.f_t;
g_f = Material.g_f;
f_c = Material.f_c;

% Previous increment
% -------------------
e_i = Material_State.e;
s_i = Material_State.s_eff;       % effective stress used by return mapping
k_Dt = Material_State.k_Dt;
k_Dc = Material_State.k_Dc;
k_i = Material_State.k;

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

% Elastic compliance reused throughout the material update.
S_e = inv(D_e);

% Elastic predictor
% ------------------
d_e = e-e_i;
s_e = s_i + D_e*d_e;

% Abaqus CDP parameters
% ----------------------
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
[F_e,~,~] = CDP_Yield(s_e,sigma_t,sigma_c,alpha,beta,gamma,tol_q);

if F_e <= tol_F
    % Elastic effective stress
    s_eff = s_e;
    D_ep  = D_e;
    k = k_i;

else
    % Return mapping: ONE yield surface + ONE plastic multiplier
    % ---------------------------------------------------------------
    s = s_e;
    k = k_i;

    converged = false;

    for it = 1:maxIter

        [F, dF_ds, ~] = CDP_Yield(s,sigma_t,sigma_c,...
                                             alpha,beta,gamma,tol_q);

        [~, dG_ds, d2G_ds2, ~] = CDP_Flow(s,sigma_t,psi,ecc,tol_q);

        % Residuals:
        %   r_sigma = De^{-1}(s-se) + (k-ki)*G_,s = 0
        %   r_F     = F(s) = 0
        dk = k-k_i;

        r1 = S_e*(s-s_e) + dk*dG_ds;
        r2 = F;
        r  = [r1;r2];

        % Dimensionally consistent convergence scaling:
        % r1 is strain-like; r2=F is stress-like.
        strain_scale = max([1.0e-12, norm(d_e), norm(S_e*s_e)]);
        stress_scale = max([1.0, f_t, f_c, norm(s_e)]);

        if norm(r1) <= tol_R*strain_scale && ...
           abs(r2) <= tol_F*stress_scale
            converged = true;
            break
        end

        % Exact-form Newton Jacobian for the chosen constitutive equations:
        %
        % dr_sigma/ds = De^{-1} + dk*d2G/ds2
        % dr_sigma/dk = dG/ds
        % dr_F/ds     = dF/ds
        % dr_F/dk     = 0  (Stage 1: no hardening evolution yet)
        J11 = S_e + dk*d2G_ds2;
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

        % Dimensionless residual-based backtracking for robustness.
        % Use the same strain/stress scaling as the convergence check so
        % strain-like r1 and stress-like r2=F are not mixed directly.
        base_norm = sqrt((norm(r1)/strain_scale)^2 + ...
                         (r2/stress_scale)^2);
        scale = 1.0;
        accepted = false;

        while scale >= 1e-6
            s_try = s + scale*ds;
            k_try = k + scale*dk_new;

            [F_try,~,~] = CDP_Yield(s_try,sigma_t,sigma_c,...
                                       alpha,beta,gamma,tol_q);
            dk_try = k_try-k_i;
            [~,dG_try,~,~] = CDP_Flow(s_try,sigma_t,psi,ecc,tol_q);
            r1_try = S_e*(s_try-s_e) + dk_try*dG_try;

            try_norm = sqrt((norm(r1_try)/strain_scale)^2 + ...
                            (F_try/stress_scale)^2);

            if try_norm < base_norm
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

    [~,dF_ds,~] = CDP_Yield(s,sigma_t,sigma_c,...
                                     alpha,beta,gamma,tol_q); %#ok<ASGLU>
    [~,dG_ds,d2G_ds2,~] = CDP_Flow(s,sigma_t,psi,ecc,tol_q);

    % H is the stress-space tangent associated with the current return
    % mapping equation.
    Hinv = S_e + dk*d2G_ds2;
    H    = inv(Hinv);

    % Yield consistency uses F_,s while plastic flow uses G_,s.
    denom = dF_ds'*H*dG_ds;

    if abs(denom) < 1e-14
        error('Damage_Plasticity_Model:SingularTangent',...
              'The unified CDP consistent tangent denominator is too small.');
    end

    D_ep = H - H*dG_ds*(1/denom)*dF_ds'*H;

end

% STAGE II-A damage update and exact branchwise nominal tangent.
% This occurs after both elastic and plastic effective-stress updates.
% For elastic steps D_ep = D_e; plastic steps use Stage-I D_ep.
% In the plastic branch s is still the effective stress at this point.
if F_e > tol_F
    s_eff = s;
end

% Abaqus-inspired stress-weighted INCREMENTAL plastic-strain histories.
% Plastic increment from the effective-stress return mapping:
%   Delta_ep = Delta_e - De^{-1} Delta_s_eff.
% Unlike Stage II-A's previous max-norm of total plastic strain, these
% histories accumulate positive increments and distinguish stress states.
% r = sum(<sigma_i_eff>_+) / sum(abs(sigma_i_eff)).
% Delta_kDt = r * <largest principal Delta_ep>_+;
% Delta_kDc = (1-r) * <-smallest principal Delta_ep>_+.
% This is an Abaqus-inspired discretization, not a claim of bitwise
% equivalence to Abaqus internal time integration.
% Elastic effective-stress steps have exactly zero plastic increment and
% zero derivative of that increment (D_ep = D_e). Therefore damage histories
% and damage gradients do not change on such steps. This avoids all history
% eigendecompositions and, for virgin material, the stress-split tangent.
if F_e <= tol_F
    d_t = 1-exp(-f_t*k_Dt/g_f);
    d_c = 0;
    % Elastic effective-stress step: damage histories are frozen, so both
    % damage derivatives are zero. Bypass the spectral split only when
    % neither tensile nor compressive damage is present.
    if d_t == 0 && d_c == 0
        s = s_eff;
        D = D_ep;
    else
        [s,D] = CDP_SplitDamageStressTangent(s_eff,D_ep,d_t,d_c,...
                                              zeros(6,1),zeros(6,1));
    end
else
    delta_ep = d_e - S_e*(s_eff-s_i);
    d_delta_ep_d_e = eye(6) - S_e*D_ep;

    % Evaluate the weight and extreme plastic strains ONCE, including their
    % gradients. The original version repeated both eigendecompositions.
    [r,grad_r_s] = CDP_TensionStressWeight(s_eff);
    [ep_max,grad_ep_max,ep_min,grad_ep_min] = ...
        CDP_ExtremePlasticStrains(delta_ep);
    pt = max(ep_max,0);
    pc = max(-ep_min,0);
    dk_t = r*pt;
    dk_c = (1-r)*pc;

    grad_pt = zeros(6,1);
    grad_pc = zeros(6,1);
    if ep_max > 0
        grad_pt = grad_ep_max;
    end
    if ep_min < 0
        grad_pc = -grad_ep_min;
    end
    grad_kDt_e = r*(d_delta_ep_d_e'*grad_pt) ...
                  + pt*(D_ep'*grad_r_s);
    grad_kDc_e = (1-r)*(d_delta_ep_d_e'*grad_pc) ...
                  - pc*(D_ep'*grad_r_s);

    k_Dt = k_Dt + dk_t;
    k_Dc = k_Dc + dk_c;
    d_t = 1-exp(-f_t*k_Dt/g_f);
    d_c = 0; % Stage II-C will introduce compression damage.
    grad_dt_e = (f_t/g_f)*(1-d_t)*grad_kDt_e;
    grad_dc_e = zeros(6,1);

    % This shortcut is valid only if damage AND its strain derivatives vanish.
    if d_t == 0 && d_c == 0 && ...
            ~any(grad_dt_e) && ~any(grad_dc_e)
        s = s_eff;
        D = D_ep;
    else
        [s,D] = CDP_SplitDamageStressTangent(s_eff,D_ep,d_t,d_c,...
                                              grad_dt_e,grad_dc_e);
    end
end

% Final material state
% ---------------------
Material_State.e     = e;
Material_State.s     = s;
Material_State.s_eff = s_eff;

% Final Stage-II-A internal variables
Material_State.k = k;

Material_State.k_Dt = k_Dt;
Material_State.k_Dc = k_Dc;
Material_State.d_t  = d_t;
Material_State.d_c  = d_c;

end


% ========================================================================
function [F,dF_ds,aux] = CDP_Yield(s,sigma_t,sigma_c,alpha,beta,gamma,tol_q)
% Abaqus CDP Lubliner-Lee-Fenves yield surface.
%
% F = 1/(1-alpha) * [ q - 3*alpha*p
%                     + beta*<sigma_max>
%                     - gamma*<-sigma_max> ] - sigma_c
%
% Only the yield gradient is required by the Stage-1 return mapping and
% consistent tangent. The unused numerical Hessian of F has been removed.

[p,q,sdev,sigma_max,Pmax] = CDP_Invariants(s,tol_q);

pos = max(sigma_max,0);
neg = max(-sigma_max,0);

F = (q - 3*alpha*p + beta*pos - gamma*neg)/(1-alpha) - sigma_c;

dp = [-1/3;-1/3;-1/3;0;0;0];

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

dsm = [Pmax(1,1);Pmax(2,2);Pmax(3,3);...
       2*Pmax(1,2);2*Pmax(2,3);2*Pmax(3,1)];

% Centered subgradient at the Macaulay-bracket corner.
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

aux.p=p; aux.q=q; aux.sdev=sdev; aux.sigma_max=sigma_max;
end


function [G,dG_ds,d2G_ds2,aux] = CDP_Flow(s,sigma_t,psi,ecc,tol_q)
% Abaqus CDP hyperbolic flow potential:
% G = sqrt((ecc*sigma_t*tan(psi))^2 + q^2) - p*tan(psi)
%
% The gradient and Hessian used by Newton and the consistent tangent are
% analytical. No finite-difference derivatives are used here.

[p,q,sdev,~,~]=CDP_Invariants(s,tol_q);

A=ecc*sigma_t*tan(psi);
R=sqrt(A^2+q^2);

G=R-p*tan(psi);

dp=[-1/3;-1/3;-1/3;0;0;0];

% q^2 = s' * Mq * s for stress Voigt order [11 22 33 12 23 31].
% The shear stresses are tensor shear components (not engineering strains).
Mq = [ 1.0 -0.5 -0.5  0 0 0;...
      -0.5  1.0 -0.5  0 0 0;...
      -0.5 -0.5  1.0  0 0 0;...
       0    0    0     3 0 0;...
       0    0    0     0 3 0;...
       0    0    0     0 0 3];

if q>tol_q
    dq=(Mq*s)/q;
    d2q=Mq/q-(dq*dq')/q;
else
    % At q=0, q itself is nondifferentiable. The hyperbolic potential G,
    % however, is smooth because A>0 for the present CDP parameters.
    dq=zeros(6,1);
    d2q=zeros(6,6);
end

dG_ds=(q/R)*dq-tan(psi)*dp;

% Exact Hessian away from q=0:
% d2G = (A^2/R^3) dq*dq' + (q/R) d2q.
% At q=0, use the finite analytical limit d2G = Mq/A.
if q>tol_q
    d2G_ds2=(A^2/R^3)*(dq*dq')+(q/R)*d2q;
else
    d2G_ds2=Mq/A;
end

% Remove roundoff-level asymmetry.
d2G_ds2=0.5*(d2G_ds2+d2G_ds2');

aux.p=p; aux.q=q; aux.A=A;
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

% Derivative/subgradient of the maximum principal stress.
%
% For a unique maximum principal stress:
%     d(sigma_max)/dS = n_max*n_max'
%
% At a repeated maximum eigenvalue the maximum-principal direction is not
% unique. Selecting only one eigenvector makes the tangent depend on the
% arbitrary eigenvector returned by eig().  Use the centered projector over
% the complete repeated maximum eigenspace instead.  This preserves the
% symmetry of axisymmetric/repeated-principal-stress states and corresponds
% to the centered derivative used by the independent tangent check.
eig_scale = max([1;abs(lam)]);
eig_tol   = 1.0e-10*eig_scale;

is_max = abs(lam-sigma_max) <= eig_tol;
m_max  = sum(is_max);

Vmax = V(:,is_max);
Pmax = (Vmax*Vmax')/m_max;

% Remove roundoff-level asymmetry.
Pmax = 0.5*(Pmax+Pmax');

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


% ========================================================================
function T = CDP_VoigtStressTensor(v)
T = [v(1),v(4),v(6);v(4),v(2),v(5);v(6),v(5),v(3)];
end

function v = CDP_TensorStressVoigt(T)
v = [T(1,1);T(2,2);T(3,3);T(1,2);T(2,3);T(3,1)];
end

function [r,grad] = CDP_TensionStressWeight(s)
% r = sum positive effective principal stresses / sum absolute stresses.
% Exact derivative on smooth stress branches, centered derivative at zero.
T=CDP_VoigtStressTensor(s);
[V,L]=eig(T);
lam=diag(L);
tol=1e-12*max(1,max(abs(lam)));
p=max(lam,0);
a=sum(p); b=sum(abs(lam));
grad=zeros(6,1);
if b<=tol
    r=0; % zero stress: no plastic loading expected
    return
end
r=a/b;
ga=zeros(3,3); gb=zeros(3,3);
for i=1:3
    if lam(i)>tol
        da=1; db=1;
    elseif lam(i)<-tol
        da=0; db=-1;
    else
        da=0.5; db=0;
    end
    ga=ga+da*(V(:,i)*V(:,i)');
    gb=gb+db*(V(:,i)*V(:,i)');
end
gr=(b*ga-a*gb)/(b*b);
grad=[gr(1,1);gr(2,2);gr(3,3);...
      2*gr(1,2);2*gr(2,3);2*gr(3,1)];
end

function [emax,gmax,emin,gmin] = CDP_ExtremePlasticStrains(dep)
% Engineering shear in strain Voigt: tensor off-diagonals are dep(4:6)/2.
T=[dep(1) dep(4)/2 dep(6)/2;...
   dep(4)/2 dep(2) dep(5)/2;...
   dep(6)/2 dep(5)/2 dep(3)];
[V,L]=eig(T);
lam=diag(L);
emax=max(lam); emin=min(lam);
tol=1e-10*max(1,max(abs(lam)));
imax=abs(lam-emax)<=tol;
imin=abs(lam-emin)<=tol;
Pmax=(V(:,imax)*V(:,imax)')/sum(imax);
Pmin=(V(:,imin)*V(:,imin)')/sum(imin);
gmax=[Pmax(1,1);Pmax(2,2);Pmax(3,3);...
      Pmax(1,2);Pmax(2,3);Pmax(3,1)];
gmin=[Pmin(1,1);Pmin(2,2);Pmin(3,3);...
      Pmin(1,2);Pmin(2,3);Pmin(3,1)];
end

function [s_nom,D_nom] = CDP_SplitDamageStressTangent(s_eff,D_ep,dt,dc,gdt,gdc)
% Exact Frechet derivative on smooth spectral branches.
% At a zero principal stress, the centered generalized derivative is used.
% Stress Voigt uses tensor shear, strain Voigt uses engineering shear.
S = CDP_VoigtStressTensor(s_eff);
[V,L] = eig(S);
lam = diag(L);
pos = max(lam,0);
neg = min(lam,0);
Sp = V*diag(pos)*V';
Sn = V*diag(neg)*V';
s_nom = CDP_TensorStressVoigt((1-dt)*Sp+(1-dc)*Sn);

% Loewner matrix for positive-part spectral function.
scale = max(1,max(abs(lam)));
Lpos = zeros(3,3);
for i=1:3
    for j=1:3
        if i==j
            if lam(i)>1e-12*scale
                Lpos(i,j)=1;
            elseif lam(i)<-1e-12*scale
                Lpos(i,j)=0;
            else
                Lpos(i,j)=0.5;
            end
        elseif abs(lam(i)-lam(j))>1e-12*scale
            Lpos(i,j)=(pos(i)-pos(j))/(lam(i)-lam(j));
        else
            if lam(i)>1e-12*scale && lam(j)>1e-12*scale
                Lpos(i,j)=1;
            elseif lam(i)<-1e-12*scale && lam(j)<-1e-12*scale
                Lpos(i,j)=0;
            else
                Lpos(i,j)=0.5;
            end
        end
    end
end

% d(sigma_nom)/d(s_eff) = (1-dc)*I + (dc-dt)*d(S_positive)/dS.
% This avoids forming dSn, redundant tensor additions and extra 3x3 products.
% Precompute the six transformed stress-Voigt basis tensors once per call.
A=(1-dc)*eye(6);
if dt ~= dc
    for j=1:6
        basis=zeros(6,1);
        basis(j)=1;
        B=V'*CDP_VoigtStressTensor(basis)*V;
        dSp=V*(Lpos.*B)*V';
        A(:,j)=A(:,j)+(dc-dt)*CDP_TensorStressVoigt(dSp);
    end
end
D_nom=A*D_ep-CDP_TensorStressVoigt(Sp)*gdt' ...
              -CDP_TensorStressVoigt(Sn)*gdc';
end
