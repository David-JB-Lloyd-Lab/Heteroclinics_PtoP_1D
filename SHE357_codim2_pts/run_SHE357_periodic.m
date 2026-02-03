
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes 2-pi periodic stationary solutions to the stationary 
% cubic quintic septic Swift-Hohenberg equation (SH357)
%         −(1 + k^2*d_zz)^2[u] + μ*u - a u^3 + bu^5 - u^7  = 0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, close all, clc;

% Spatial coordinates: z direction
nz = 40; Lz = pi; hz = 2*pi/nz;  z = hz*(1:nz); z = Lz*(z-pi)/pi;
% Fourier psuedo-spectral differentiation matrices
column = [0 .5*(-1).^(1:nz-1).*cot((1:nz-1)*hz/2)]';
Dz  = toeplitz(column,column([1 nz:-1:2]));
D2z = toeplitz([-pi^2/(3*hz^2)-1/6 .5*(-1).^(2:nz)./sin(hz*(1:nz-1)/2).^2]);
D4z = D2z^2;
wz  = 2*pi*ones(1,nz)/nz; % integration weights for trapzoid rule - mean
z   = z';
Iz = speye(nz);

mesh_params.nz  = nz;  mesh_params.Lz  = Lz;  mesh_params.z   = z;    
mesh_params.Iz  = Iz;  mesh_params.wz  = wz;  
mesh_params.Dz  = Dz;  mesh_params.D2z = D2z; mesh_params.D4z = D4z;

 % Swift-Hohenberg equation parameters p = [mu,a,b,sigma,k,gamma];
p(1) = 0.3; %mu 
p(2) = 1.5; %a 
p(3) = 1.5; %b 
p(4) = 2; %sigma 
p(5) = 0.8; %k 
p(6) = 1; %gamma  

% initial guess
w0 = 1*cos(z(:));
mesh_params.w0 = w0; % reference profiles
mesh_params.w0z= Dz*w0;
u0 = [w0; 0];

%% converge to an initial periodic solution
disp('Converging periodic orbit up');
my_rhs = @(u) SHE357_1D(u,p,mesh_params);
options = optimset('Jacobian','on','Display','iter','MaxIter',500,'TolFun',1e-12);
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

%% Continue in the bifurcation parameter. Continuation is from: 
% D. Avitabile, "Numerical Computation of Coherent Structures in
% Spatially-Extended Systems", Zenodo (2020).
% DOI: 10.5281/zenodo.3821169
problemHandle            = @(u,p)SHE357_1D(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_SHE357_periodic(step,u,p,mesh_params);
computeEigenvaluesHandle = @(u,p) ComputeSpectrum(u,p,mesh_params);
plotSpetcrumHandle       = @(d,p,parentHandle) PlotSpectrum(d,p,parentHandle); 
stepperPars.iContPar      = 5; % parameter to continue 
stepperPars.s0            = 0.01; %or -0.01
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 0.01;
stepperPars.pMin          = sqrt(1 - sqrt(p(1)))  + 1e-3;
stepperPars.pMax          =  1;
stepperPars.maxSteps      = 200;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'Jacobian','on',...
                                     'MaxIter',30,'TolFun',1e-10);
stepperPars.optNonlinIter = 20;
stepperPars.dataFolder    = 'periodic_solutions'; 
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = []; plotSpetcrumHandle; 
stepperPars.PlotBranchVariableId = 3; % plot this index of sol measures
branch = SecantContinuation(problemHandle,u_out,p,stepperPars);

%save mesh parameters. If the continuations is stopped prematurely, this
%may have to be done manually. 
if ~exist(stepperPars.dataFolder, 'dir')
    mkdir(stepperPars.dataFolder)
end
save(fullfile(stepperPars.dataFolder, 'mesh_params'), 'mesh_params');
pause(0.01)