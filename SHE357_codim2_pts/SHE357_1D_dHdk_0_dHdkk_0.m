function [F] = SHE357_1D_dHdk_0_dHdkk_0(uin,p,mesh_params,~)

% Rename parameters
mu = p(1);
a  = p(2);
b  = p(3);
sigma = p(4);
kx = p(5);

% Auxiliary variables
n = mesh_params.nz;
u = uin(1:n);
uk= uin(1+n:2*n);
ukk = uin(2*n+1:3*n);
c1 = uin(3*n+1);
c2 = uin(3*n+2);
c3 = uin(3*n+3);
b = uin(3*n+4);
kx = uin(3*n + 5);

D4z = mesh_params.D4z; % differentiation matrices
D2z = mesh_params.D2z;
Dz  = mesh_params.Dz;

uxxxx= (kx)^4*D4z*u;
uxx  = kx^2*D2z*u;
ux   = kx*Dz*u;

u4xk = (kx)^4*D4z*uk;
u2xk = kx^2*D2z*uk;
uxk  = kx*Dz*uk;

u4xkk = (kx)^4*D4z*ukk;
u2xkk = kx^2*D2z*ukk;
uxkk  = kx*Dz*ukk;

% SHE
F1 = -uxxxx - sigma*uxx - u + mu*u - a*u.^3 + b*u.^5 - u.^7 + c1*ux;

% Derivative of periodic orbit with respect to k
F2  = -(4/kx)*uxxxx - u4xk  ...
    - (2/kx)*sigma*uxx - sigma*u2xk ...
    - uk + mu*uk - 3*a*u.^2.*uk + 5*b*u.^4.*uk - 7*u.^6.*uk + c2*uxk; 

% Second derivative of periodic orbit with respect to k
F3 = -(12/kx^2)*uxxxx - (8/kx)*u4xk - u4xkk ...
    - (2/kx^2)*sigma*uxx - sigma*(4/kx)*u2xk - sigma*u2xkk ...
    - ukk + mu*ukk ...
    - 6*a*u.*uk.^2 - 3*a*u.^2.*ukk + 20*b*u.^3.*uk.^2 + 5*b*u.^4.*ukk ...
    - 42*u.^5.*uk.^2 - 7*u.^6.*ukk - c3*uxkk;

loc = mesh_params.w0z;
F4 =  mesh_params.wz*(loc.*(u  - mesh_params.w0)); % phase conditions
F5 =  mesh_params.wz*(loc.*(uk  - mesh_params.w0));
F6 =  mesh_params.wz*(loc.*(ukk  - mesh_params.w0));

H = -0.5*kx^4*(D2z*u).^2 + sigma*kx^2*(Dz*u).^2/2 + (1-mu)*u.^2/2  ...
    + 1/4*a*u.^4 - b*u.^6/6 + u.^8/8 + kx^4*(Dz*u).*(Dz*D2z*u); % Hamiltonian

% dH/dk
dHdk = -2*kx^3*(D2z*u).^2 -kx^4*(D2z*u).*(D2z*uk) ...
    + sigma*kx*(Dz*u).^2 + sigma*kx^2*(Dz*u).*(Dz*uk) ...
    + (1-mu)*u.*uk + a*u.^3.*uk - b*u.^5.*uk + u.^7.*uk ...
    + 4*kx^3*(Dz*u).*(Dz*D2z*u) + kx^4*(Dz*uk).*(Dz*D2z*u) ...
    + kx^4*(Dz*u).*(Dz*D2z*uk);

k = kx;

% d^2H/dk^2
dHdkk = ...
    -6*k^2.*(D2z*u).^2 ...
    -8*k^3.*(D2z*u).*(D2z*uk) ...
    - k^4.*((D2z*uk).^2 + (D2z*u).*(D2z*ukk)) ...
    + sigma.*(Dz*u).^2 + 4*sigma.*k.*(Dz*u).*(Dz*uk) + sigma.*k^2.*((Dz*uk).^2 + (Dz*u).*(Dz*ukk)) ...
    + (1-mu).*( uk.^2 + u.*ukk ) ...
    + a.*( 3*u.^2.*uk.^2 + u.^3.*ukk ) ...
    - b.*( 5*u.^4.*uk.^2 + u.^5.*ukk ) ...
    + ( 7*u.^6.*uk.^2 + u.^7.*ukk ) ...
    + 12*k^2.*(Dz*u).*(Dz*D2z*u) ...
    + 8*k^3.*(Dz*uk).*(Dz*D2z*u) ...
    + 8*k^3.*(Dz*u).*(Dz*D2z*uk) ...
    + k^4.*( (Dz*ukk).*(Dz*D2z*u) + 2*(Dz*uk).*(Dz*D2z*uk) + (Dz*u).*(Dz*D2z*ukk) );

F7 = mean(dHdk); % set dH/dk =0
F8 = mean(dHdkk); % set d^2H/dk^2 = 0

F = [F1;F2;F3;F4;F5;F6;F7;F8];

end
