function [uu_1,Duu_1,uu_2,Duu_2,uu] = get_rolls(p,k,phi,mesh_params,ur_0)
  
  % construct differential operators operators
  Dy_scale = mesh_params.Dy*k;
  Dyy_scale= mesh_params.D2y*k^2;
  
  % Solve SH on a periodic domain 
  ny = mesh_params.ny; hy = 2*pi/ny;  y = hy*(1:ny)';

  % find 1D rolls
  RD_rhs_1D = @(u) rhs_rolls_phaseCond(u,p,Dyy_scale,Dy_scale,y,mesh_params);
  options = optimset('Jacobian','on','Display','off');
  [ur] = fsolve(RD_rhs_1D,[ur_0;0],options);
  
  if abs(ur(end))>1e-3
      disp(['Warning phase condition for rolls is non-zero =',num2str(ur(end))]);
  end
  
  uu = ur(1:end-1);
  
  u1 = ur(1:ny); % remove the wave speed from ur
  u2 = ur(1+ny:2*ny);
  
  Du1 = Dy_scale*u1;
  Du2 = Dy_scale*u2;
  
  % Interpolating on the eta mesh
  eta = k*mesh_params.x + phi;  % find skewed rolls and interpolate

  uu_1 = zeros(size(eta));
  Duu_1= zeros(size(eta));
  uu_2 = zeros(size(eta));
  Duu_2= zeros(size(eta));
  periodic_sinc = ones(size(eta)); 

for i = 1:ny
    xx = eta-y(i);
    periodic_sinc = diric(xx,ny).*cos(xx/2);
    uu_1 = uu_1 + u1(i)*periodic_sinc;
    Duu_1= Duu_1+Du1(i)*periodic_sinc;
    uu_2 = uu_2 + u2(i)*periodic_sinc;
    Duu_2= Duu_2+Du2(i)*periodic_sinc;
end

