function F = SolutionMeasures_SHE357_periodic(step,u_out,p,mesh_params)
  
mu = p(1);
a = p(2);
b  = p(3);
sigma = p(4);
k  = p(5);
g = p(6); 

  % Auxiliary variables
  n = mesh_params.nz;
  u = u_out(1:n);
  
  w = mesh_params.wz;
  D2z = mesh_params.D2z;
  Dz = mesh_params.Dz;

  I = speye(n);


  E = w*((k^2.*D2z*u).^2 - sigma*(k*Dz*u).^2/2 - 0.5*mu*u.^2 + a*u.^4/4 - b*u.^6/6 + u.^8/8);
  
   H = -0.5*k^4*(D2z*u).^2 + sigma*k^2*(Dz*u).^2/2 + (1-mu)*u.^2/2  ...
      + 1/4*a*u.^4 - b*u.^6/6 + u.^8/8 + k^4*(Dz*u).*(Dz*D2z*u);
  
  A = mean(2*k*(Dz*u).^2 - 2*k^3*(D2z*u).^2); %Action

  F = [mean(u.^2), E, mean(H) , A]; 

end



function [F] = SHE357_1D_get_Hk(uk,uin,p,mesh_params,~)

% Rename parameters
mu = p(1);
a  = p(2);
b  = p(3);
sigma = p(4);

% Auxiliary variables
n = mesh_params.nz;
u = uin(1:n);
kx = uin(end-1);
c = uin(end); 

D4z = mesh_params.D4z;
D2z = mesh_params.D2z;
Dz  = mesh_params.Dz;

uxxxx= (kx)^4*D4z*u;
uxx  = kx^2*D2z*u;
ux   = kx*Dz*u;

u4xk = (kx)^4*D4z*uk;
u2xk = kx^2*D2z*uk;
uxk  = kx*Dz*uk;

loc = mesh_params.w0z; 
F1 = -(4/kx)*uxxxx - u4xk  ...
    - (2/kx)*sigma*uxx - sigma*u2xk ...
    - uk + mu*uk - 3*a*u.^2.*uk + 5*b*u.^4.*uk - 7*u.^6.*uk + c*uxk; 
F2 =  mesh_params.wz*(loc.*(uk  - mesh_params.w0));


F = [F1; F2]; 
end
