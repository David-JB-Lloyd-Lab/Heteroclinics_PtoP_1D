function [F,J] = SHE357_1D(uin,p,mesh_params)

  % Rename parameters
mu = p(1);
a = p(2);
b  = p(3);
sigma = p(4);
kx = p(5);
gamma = p(6); 
  
  % Auxiliary variables
  n = mesh_params.nz;
  u = uin(1:n);
  c = uin(n+1);

  uxxxx= (kx)^4*mesh_params.D4z*u;
  uxx  = kx^2*mesh_params.D2z*u;
  ux   = kx*mesh_params.Dz*u;
  
  % SHE equation
  F = -uxxxx - sigma*uxx - u + mu*u - a*u.^3 + b*u.^5 - gamma*u.^7 + c*ux;
  
  loc = mesh_params.w0z;
  F(n+1) =  mesh_params.wz*(loc.*(u  - mesh_params.w0)); % phase condition
  
  % Jacobian
  if nargout > 1
     J = sparse(n,n);
     J = -kx^4*mesh_params.D4z - sigma*kx^2*mesh_params.D2z ...
         + spdiags((-1+mu)*ones(n,1) - 3*a*u.^2 + 5*b*u.^4 - gamma*7*u.^6,0,n,n) ...
         + c*kx*mesh_params.Dz;
%      
     J(n+1,1:n)  =mesh_params.wz*spdiags(loc,0,n,n);

     epsiF = 1e-8;
     dF = SHE357_1D([uin(1:n); uin(n+1) + epsiF],p,mesh_params);
     J(:,n+1) = (dF - F)/epsiF;

  end
      

end
