%computation of stationary PtoP solutions for the ac-boussinesq equation 
% a u_xx + u + h(u + \alpha u^2)          = A 
% b h_xx + h + 0.5(h^2) + (1/3)\alpha u^3 = B
%using the core-farfield decomposition method. 

%this script computes periodic orbits in the far-field u_\pm(x; k_\pm,\phi_\pm) and
%the core of the solution. The contination is designed so that the front
%approaches the codimension-2 point where the wavenumbers of the farfield
%orbits coalesce

clear all, close all, clc;
global ur0_left ur0_right; % make global far-field periodic orbit

% a-c Boussinesq parameters
a = 0.2;
b = 0; %coefficient in the Boussinesq models proposed by Bona, Chen, Saut 
c = -0.3;
A = -2;
B = 2;

kxm = 1.967; % initial guess for selected wavenumber at x = -infty
kxp = 2.2797; % initial guess for selected wavenumber at x = +infty
d = -60;    % cutoff point for partion function
m = 2;      % stepness of partion funtion
phi = 0;    % asymptotic phase
alpha = 0.1935; 

p(1)  = a;
p(2)  = b;
p(3)  = c;
p(4)  = A;
p(5)  = B;
p(6)  = 0;
p(7)  = kxm;
p(8)  = kxp; 
p(9)  = alpha;
p(10) = m;  
p(11) = d;
p(12) = phi;

%% define nonlinearity and Jacobian and store in mesh_params structure
mesh_params.ff = @(u,p) 0.5*u.^2 + p(9)/3*u.^3 ;
mesh_params.ffp = @(u,p) u + p(9)*u.^2 ;
mesh_params.ffpp = @(u,p) 1 + 2*p(9)*u ;

mesh_params.F1 = @(u,h,p) u + h.*mesh_params.ffp(u,p) ...
    - p(4)*ones(size(h));
mesh_params.F2 = @(u,h,p) h + mesh_params.ff(u,p) ...
    - p(5)*ones(size(u));

mesh_params.DF11 = @(u,h,p) ones(size(u)) + h.*mesh_params.ffpp(u,p); % Jacobian
mesh_params.DF12 = @(u,h,p) mesh_params.ffp(u,p) - p(6)*ones(size(h));
mesh_params.DF21 = @(u,h,p) -p(6) + mesh_params.ffp(u,p);
mesh_params.DF22 = @(u,h,p) ones(size(h));

%% Setup spatial mesh
% Spatial coordinates: x direction
nx = 2^12; Lx = 150; x = linspace(-Lx/2,Lx/2,nx)'; hx = x(2)-x(1);
[~,D2x,Dx,wx,D4x] = Compute_1D_5ptLaplacian_finite_difference(nx,Lx);

% Spatial discretisation in y
ny = 32; Ly = pi; hy = 2*pi/ny;  y = hy*(1:ny)'; y = Ly*(y-pi)/pi;
[~,D2y,Dy,wy] = Compute_1D_Laplacian_fourier(ny,Ly);

% store mesh and differentiation matrices in mesh_params structure
mesh_params.Dx = Dx; mesh_params.D2x = D2x; mesh_params.D4x = D4x; mesh_params.wx = wx;
mesh_params.nx = nx; mesh_params.Lx  = Lx;  mesh_params.x   = x;
mesh_params.y = y; mesh_params.Ly = Ly; mesh_params.ny = ny;
mesh_params.Dy = Dy; mesh_params.D2y = D2y; mesh_params.wy = wy;
ii = find(abs(x) > Lx/2-4*hx); mesh_params.iiend = ii;

%% setup initial condition for far-field single periodic orbit
u0 = -2 + cos(y);
h0 = 1.5 + 0.5*cos(y);
ur0_left = [ u0; h0];
mesh_params.uT1 = ur0_left(1:ny); mesh_params.uT2 = ur0_left(1+ny:end);
[uuL_1,~,uuL_2,~,ur0_left] = get_rolls(p,kxm,0,mesh_params,ur0_left); 

u0 = -1.43 + cos(y);
h0 =  1 + 0.5*cos(y);
ur0_right = [u0; h0];
mesh_params.uT1 = ur0_right(1:ny); mesh_params.uT2 = ur0_right(1+ny:end);
[uuR_1,~,uuR_2,~,ur0_right] = get_rolls(p,kxp,phi,mesh_params,ur0_right); 
% get/converge far-field periodic orbit


%% construct initial guess for front
% partion functions
chi_p = 1/2 + 1/2*tanh(m*(x-d)); chi_m = chi_p(end:-1:1);
c1p = 1/2*(1+tanh(1*x)); c1m = c1p(end:-1:1);
w0 = [uuL_1.*chi_p.*c1m + uuR_1.*chi_m.*c1p;
    uuL_2.*c1m.*chi_p + uuR_2.*chi_m.*c1p];

u0 = [w0; kxm; kxp]; %initial guesss 

% converge initial guess
my_rhs = @(u) Bouss_Front_fix_phase(u,p,mesh_params,1);
options = optimset('Jacobian','on','Display','off','MaxIter',5e5 ...
   ,'Algorithm','Levenberg-Marquardt','TolX',100*eps,'TolFun',1e-10);
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

%% continue
problemHandle            = @(u,p) Bouss_Front_fix_phase(u,p,mesh_params,1);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_BoussinesqFront(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
stepperPars.iContPar      = 9; % cont parameter index
stepperPars.s0            = 0.001; % initial paramter step (-0.1 to step in the opposite direction)
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 1;
stepperPars.pMin          = -10;
stepperPars.pMax          = 10;
stepperPars.maxSteps      = 1000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
options = optimset('Jacobian','on','Display','off','MaxIter',100); % fsolve options
options.border = nx;
stepperPars.fsolveOptions = options;
stepperPars.optNonlinIter = 50;
stepperPars.dataFolder    = 'test';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;      
stepperPars.PlotBranchVariableId = 2;
branch = SecantContinuation(problemHandle,u_out,p,stepperPars);