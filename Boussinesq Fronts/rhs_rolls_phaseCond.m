function [F,J] = rhs_RD_phaseCond(uu,p,Dyy,Dy,y,mesh_params)

  % Rename parameters
  n  = length(y);
  a  = p(1);
  b  = p(2);
  c  = p(3);

  u = uu(1:n);
  h = uu(1+n:2*n);
  v = uu(2*n+1);

  % Right-hand side
  uT1 = mesh_params.uT1;
  uT2 = mesh_params.uT2;
  F = zeros(2*n+3,1);

  % RD system
  F1= a*(Dyy*u) - v*Dy*u + mesh_params.F1(u,h,p);
  F2= c*(Dyy*h) - v*Dy*h + mesh_params.F2(u,h,p);
  F3 = mean((Dy*uT1).*(u - uT1) + (Dy*uT2) .* ( h - uT2)) ; % phase condition
  
  F = [F1; F2; F3];

  % Jacobian
  if nargout > 1
    J = sparse(2*n+3);
    e = ones(n,1);
    J11 = spdiags(mesh_params.DF11(u,h,p),0,n,n) + a*Dyy - v*Dy;
    J12 = spdiags(mesh_params.DF12(u,h,p),0,n,n);
    J13 = -Dy*u;
    J21 = spdiags(mesh_params.DF21(u,h,p),0,n,n);
    J22 = spdiags(mesh_params.DF22(u,h,p),0,n,n) + c*Dyy - v*Dy;
    J23 = -Dy*h;
    J31 = (Dy*uT1)'/n;
    J32 = (Dy*uT2)'/n;
    J33 = sparse(0);
     
    J = [J11 J12 J13
         J21 J22 J23
         J31 J32 J33];

  end

end
