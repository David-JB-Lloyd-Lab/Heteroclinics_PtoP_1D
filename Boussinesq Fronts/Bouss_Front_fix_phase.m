function [F,J] = Bouss_Front_fix_phase(u,p,mesh_params,jac)

global ur0_left ur0_right;

% Rename parameters
a  = p(1);
b  = p(2);
c  = p(3);
m  = p(10);
d  = p(11);
phi = p(12); 

% Auxiliary variables
n  = mesh_params.nx;
w1 = u(1:n);
w2 = u(1+n:2*n);
kxm= u(2*n+1);
kxp = u(2*n+2);

%
chi_p = 1/2 + 1/2*tanh(m*(mesh_params.x+d));
chi_m = chi_p(end:-1:1);

% find far-field periodic orbits 
if jac ~= 1
    [uuL_1,DuuL_1,uuL_2,DuuL_2,ur0_left] = get_rolls(p,kxm,0,mesh_params,ur0_left);
    [uuR_1,DuuR_1,uuR_2,DuuR_2,ur0_right] = get_rolls(p,kxp,phi,mesh_params,ur0_right);
else
    [uuL_1,DuuL_1,uuL_2,DuuL_2,~] = get_rolls(p,kxm,0,mesh_params,ur0_left);
    [uuR_1,DuuR_1,uuR_2,DuuR_2,~] = get_rolls(p,kxp,phi,mesh_params,ur0_right);
end

% Right-hand side
L1 = a*mesh_params.D2x;
L2 = c*mesh_params.D2x;

CHI_p = spdiags(chi_p,0,n,n);
CHI_m = spdiags(chi_m,0,n,n);


F = zeros(length(u),1);
%build odes for w
remFun1 = uuL_1.*chi_m + uuR_1.*chi_p + w1; 
remFun2 = uuL_2.*chi_m + uuR_2.*chi_p + w2;  
F1= L1*(remFun1) + mesh_params.F1(remFun1, remFun2,p) ...
    - chi_m.*(L1*uuL_1 + mesh_params.F1(uuL_1,uuL_2,p)) ... 
    - chi_p.*(L1*uuR_1 + mesh_params.F1(uuR_1,uuR_2,p));

F2= L2*(remFun2) + mesh_params.F2(remFun1, remFun2,p) ...
    - chi_m.*(L2*uuL_2 + mesh_params.F2(uuL_1,uuL_2,p))...
    - chi_p.*(L2*uuR_2 + mesh_params.F2(uuR_1,uuR_2,p));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase condition on left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iBm = find( mesh_params.x <= -mesh_params.Lx/2 + 2*pi/abs(kxm) ); 
% find indices close to x = Lx - 2*pi/kx : Lx and y = 0..Ly
x = linspace(-mesh_params.Lx/2,mesh_params.Lx/2,mesh_params.nx);  hx= x(2) - x(1);
iim= find(x <= -mesh_params.Lx/2 + 2*pi/abs(kxm)); nnxm = length(iim);

% phase condition 2: find du_r
uuu_rm_prime1 = DuuL_1;
uuu_rm_prime2 = DuuL_2;

% define uuu_r_prime for integral
uuu_rm_prime1 = uuu_rm_prime1(iBm);
uuu_rm_prime2 = uuu_rm_prime2(iBm);

wxm = [1, 2*ones(1,nnxm-2)+2*mod([1:nnxm-2],2),1]*hx/3;   % Simpson weights for intergration int = w*u

wBm1 = w1(iBm); % find u on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
wBm2 = w2(iBm); % find h on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
F3 = wxm*(uuu_rm_prime1.*wBm1 + uuu_rm_prime2.*wBm2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase condition on right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iBp = find( mesh_params.x(:) >= mesh_params.Lx/2 - 2*pi/abs(kxp) ); % find indices close to x = Lx - 2*pi/kx : Lx and y = 0..Ly
x = linspace(-mesh_params.Lx/2,mesh_params.Lx/2,mesh_params.nx); 
hx= x(2) - x(1);
iip= find(x >= mesh_params.Lx/2 - 2*pi/abs(kxp)); nnxp = length(iip);

% phase condition 2: find du_r
uuu_rp_prime1 = DuuR_1;
uuu_rp_prime2 = DuuR_2;

% define uuu_r_prime for integral
uuu_rp_prime1 = uuu_rp_prime1(iBp);
uuu_rp_prime2 = uuu_rp_prime2(iBp);
wxp = [1, 2*ones(1,nnxp-2)+2*mod([1:nnxp-2],2),1]*hx/3;   % Simpson weights for intergration int = w*u
wBp1 = w1(iBp); % find u on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
wBp2 = w2(iBp); % find h on the domain x = Lx - 2*pi/kx : Lx and y = 0..Ly
F4 = wxp*(uuu_rp_prime1.*wBp1 + uuu_rp_prime2.*wBp2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RHS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = [F1; F2; F3; F4]; %RHS for w

% Jacobian
if nargout > 1
    J = sparse(length(u),length(u));

    J11 = spdiags(mesh_params.DF11(remFun1, remFun2,p),0,n,n) + L1;
    J12 = spdiags(mesh_params.DF12(remFun1, remFun2,p),0,n,n);
    J21 = spdiags(mesh_params.DF21(remFun1, remFun2,p),0,n,n);
    J22 = spdiags(mesh_params.DF22(remFun1, remFun2,p),0,n,n) + L2;

    J = [J11 J12
        J21 J22];

    J(2*n+1,iBm)  = wxm*spdiags(uuu_rm_prime1,0,nnxm,nnxm);
    J(2*n+2,iBp)  = wxp*spdiags(uuu_rp_prime1,0,nnxp,nnxp);
    J(2*n+1,iBm+n)  = wxm*spdiags(uuu_rm_prime2,0,nnxm,nnxm);
    J(2*n+2,iBp+n)  = wxp*spdiags(uuu_rp_prime2,0,nnxp,nnxp);

    %Terms in Jacobain from phase condition
    epsiF = 1e-6;
    dF = Bouss_Front_fix_phase([u(1:2*n); u(2*n+1) + epsiF; u(2*n+2)],p,mesh_params,1);
    dF2 = Bouss_Front_fix_phase([u(1:2*n); u(2*n+1) - epsiF; u(2*n+2)],p,mesh_params,1);
    J(:,2*n+1) = 0.5*(dF - dF2)/epsiF;

    dF = Bouss_Front_fix_phase([u(1:2*n); u(2*n+1); u(2*n+2) + epsiF],p,mesh_params,1);
    dF2 = Bouss_Front_fix_phase([u(1:2*n); u(2*n+1); u(2*n+2) - epsiF],p,mesh_params,1);

    J(:,2*n+2) = 0.5*(dF - dF2)/epsiF;

end

end

