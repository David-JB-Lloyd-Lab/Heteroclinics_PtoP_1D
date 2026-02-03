% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function [V,LAMBDA] = ComputeSpectrum_hkk_0(uin,p,mesh_params)

  %% Compute linear operators
  mu = p(1);
  a = p(2);
  % b  = p(3);
  sigma = p(4);
  % kx = p(5);
  gamma = p(6); 
  
  % Auxiliary variables
  n = mesh_params.nz;
  u = uin(1:n);
  c = uin(n+1);
  b = uin(end-1); 
  kx = uin(end); 

  uxxxx= (kx)^4*mesh_params.D4z*u;
  uxx  = kx^2*mesh_params.D2z*u;
  ux   = kx*mesh_params.Dz*u;

  J0 = -kx^4*mesh_params.D4z - sigma*kx^2*mesh_params.D2z ...
      + spdiags((-1+mu)*ones(n,1) - 3*a*u.^2 + 5*b*u.^4 - gamma*7*u.^6,0,n,n); 
  J4 = -speye(n);
  J3 = -4*kx*mesh_params.Dz;
  J2 = -sigma*speye(n)-6*kx.^2.*mesh_params.D2z;
  J1 = -4*kx.^3.*mesh_params.Dz*mesh_params.D2z - 2*kx*sigma*mesh_params.Dz;

  %% Call direct eigenvalue solver
  [V,LAMBDA] = polyeig(J0,J1,J2,J3,J4);

  [~,ii] = sort(abs(LAMBDA));
  LAMBDA = LAMBDA(ii);
  LAMBDA = exp(2*pi*LAMBDA(1:4));

end
