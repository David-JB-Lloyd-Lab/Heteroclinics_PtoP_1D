function [F] = SHE357_1D_dHdk_0(uin,p,mesh_params,~)

  % Rename parameters
mu = p(1);
a  = p(2);
b  = p(3);
% b = 3.5 + 0.4*(a-3); 
sigma = p(4);
gamma = p(6); 
  
  % Auxiliary variables
  n = mesh_params.nz;
  u = uin(1:n);
  uk= uin(1+n:2*n);
  c1 = uin(2*n+1);
  c2 = uin(2*n+2);
  kx = uin(2*n+3);

  D4z = mesh_params.D4z;
  D2z = mesh_params.D2z;
  Dz  = mesh_params.Dz;

  uxxxx= (kx)^4*D4z*u;
  uxx  = kx^2*D2z*u;
  ux   = kx*Dz*u;

  u4xk = (kx)^4*D4z*uk;
  u2xk = kx^2*D2z*uk;
  uxk  = kx*Dz*uk;
  % SHE equation
  F1 = -uxxxx - sigma*uxx - u + mu*u - a*u.^3 + b*u.^5 - gamma*u.^7 + c1*ux;

  % derivative of periodic orbit with respect to k
  F2 = -(4/kx)*uxxxx - u4xk  ...
      - (2/kx)*sigma*uxx - sigma*u2xk ...
      - uk + mu*uk - 3*a*u.^2.*uk + 5*b*u.^4.*uk - gamma*7*u.^6.*uk + c2*kx*Dz*uk;

  loc = mesh_params.w0z;
  F3 =  mesh_params.wz*(loc.*(u  - mesh_params.w0)); % phase conditions
  F4 =  mesh_params.wz*(loc.*(uk  - mesh_params.w0)); % phase conditions

  % Hamiltonian
  H = -0.5*kx^4*(D2z*u).^2 + sigma*kx^2*(Dz*u).^2/2 + (1-mu)*u.^2/2  ...
      + 1/4*a*u.^4 - b*u.^6/6 + gamma*u.^8/8 + kx^4*(Dz*u).*(Dz*D2z*u);

  % dH/dk
  dHdk = -2*kx^3*(D2z*u).^2 -kx^4*(D2z*u).*(D2z*uk) ...
       + sigma*kx*(Dz*u).^2 + sigma*kx^2*(Dz*u).*(Dz*uk) ...
       + (1-mu)*u.*uk + a*u.^3.*uk - b*u.^5.*uk + gamma*u.^7.*uk ...
       + 4*kx^3*(Dz*u).*(Dz*D2z*u) + kx^4*(Dz*uk).*(Dz*D2z*u) ...
       + kx^4*(Dz*u).*(Dz*D2z*uk);

  % derivative of Action - dA/dk
  Ak = mean(sigma*(Dz*u).^2 + 2*sigma*kx.*(Dz*u).*(Dz*uk) ...
      - 6*kx^2*(D2z*u).^2 - 4*kx.^3*(D2z*u).*(D2z*uk));

  F5 = Ak;
  
  F = [F1;F2;F3;F4;F5];

end
