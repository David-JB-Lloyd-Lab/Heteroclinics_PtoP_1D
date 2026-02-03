function plotHandle = PlotSolution(u,p,parentHandle,mesh_params)

global ur0_left;

   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[2/4*scrsz(3) scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
   else
     plotHandle = parentHandle;
   end 

   figure(parentHandle); %hold on;
   
    A  = p(1);
    R  = p(2);
    B  = p(3);
    D1 = p(4);
    D2 = (A/R)^2;
    m  = p(5);
    d  = p(6);
    phi= p(7);
  
  % Auxiliary variables
  n = mesh_params.nx;
  
  chi_p = 1/2 + 1/2*tanh(m*(mesh_params.x+d));
  chi_m = chi_p(end:-1:1);
  kxm= u(2*n+1);
  
  [uu_1,Duu_1,uu_2,Duu_2,~       ] = get_rolls(p,kxm,phi,mesh_params,ur0_left);
  
  w1     = u(1:n);
  w2     = u(1+n:2*n);
  u1 = 0*uu_1.*chi_m + w1;
  u2 = 0*uu_2.*chi_m + w2;
 
plot(mesh_params.x,u1,'b',mesh_params.x,u2,'r'); drawnow;hold off;

   % print -dtiff state.tiff

end
