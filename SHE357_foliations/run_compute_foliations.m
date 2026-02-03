clear all; close all; clc; 
 %this script computes leafs of trajectories in the stable (unstable)
 %directions for the cubic-quadratic SH equation 
 %the stable (unstable) directions are found by perturbing in the direction
 %of the stable (unstable) eigenvector of the Monodromy matrix. The phase
 %is vareied across the periodic orbit to generate the foliation from the
 %compution of individual leafs. 

N = 2^6; %number of points. 
p(1) = -1; %speed (coefficient of u term). 
p(2) =  1.38; %wavenumber; 
p(4) = 1; %alpha (cofficient of third-order dispersion term)
p(5) = -1; %gamma (coefficient of cubic NL) 
p(6) = 0; %coeff of the quadratic NL term 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%compute periodic solution and its unstable manifold 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 2*pi/p(2)/N*[-N/2:N/2-1]'; 
u = 1 +  0.5*cos(p(2)*x); 
fsolve_opts = optimset('display','off','TolX',10*eps,'TolFun',1e-12, ... 
    'MaxIter',1e3,'MaxFunEvals',1e3); 
[U,H]  = get_per(u,p,fsolve_opts); 

%compute the monodromy matrix 
phi = get_monodromy_cubic(U,p); 
fprintf('real part of eigenvalues of monodromy matrix: %0.3f %0.3f %0.3f %0.3f \n'...
    ,eig(reshape(deval(phi,2*pi/p(2)),[4,4]))'); 

%compute multiple trajectories, and their intersection with the Poincare
nPts = 200; 
pertSize = -1e-8; 

%manifold sign (+ for unstable manifolds, - for stable manifold)
manifSign = 1; 
tf = manifSign*20*pi/p(2); %decide how far you want to integrate (if the sol grows too fast, decrease)

%define the Poincare seciton u(1) = u, u(2) = u', u(3) = u'', u(4) = u''' 
pSecfun = @(u)(u(1,:));

phase = 0; 
for phase = linspace(0,2*pi/p(2),100)
[x,y] = compute_leaf(U,p,phi,tf,pertSize,phase,nPts); 
plot3(y(:,1),y(:,2),y(:,3),'b'); hold on 
drawnow(); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%compute periodic solution to the right and its stable manifold 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 2*pi/p(2)/N*[-N/2:N/2-1]'; 
u = -(1 +  0.5*cos(p(2)*x)); 
fsolve_opts = optimset('display','off','TolX',10*eps,'TolFun',1e-12, ... 
    'MaxIter',1e3,'MaxFunEvals',1e3); 
[U,H]  = get_per(u,p,fsolve_opts); 

%compute the monodromy matrix 
phi = get_monodromy_cubic(U,p); 
fprintf('real part of eigenvalues of monodromy matrix: %0.3f %0.3f %0.3f %0.3f \n'...
    ,eig(reshape(deval(phi,2*pi/p(2)),[4,4]))'); 

%compute multiple trajectories, and their intersection with the Poincare
nPts = 200; 
pertSize = 1e-8; 

%manifold sign (+ for unstable manifolds, - for stable manifold)
manifSign = -1; 
tf = manifSign*20*pi/p(2); %decide how far you want to integrate (if the sol grows too fast, decrease)

%define the Poincare seciton u(1) = u, u(2) = u', u(3) = u'', u(4) = u''' 
pSecfun = @(u)(u(1,:));

phase = 0; 
for phase = linspace(0,2*pi/p(2),100)
[x,y] = compute_leaf(U,p,phi,tf,pertSize,phase,nPts); 
plot3(y(:,1),y(:,2),y(:,3),'r'); hold on 
drawnow(); 
end