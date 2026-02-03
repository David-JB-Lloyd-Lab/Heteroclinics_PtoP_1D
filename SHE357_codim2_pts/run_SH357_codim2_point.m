%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and continues codimension two points for periodic solutions of
% the SHE357 equation 
%         −(1 + k^2*d_zz)^2[u] + μ*u - a u^3 + bu^5 - u^7  = 0 
% with the additional constraint H_k = H_kk = 0, where H is the spatial
% Hamiltonian.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load initial guess and mesh 
load('periodic_solutions_hk_0/mesh_params.mat'); 
load("periodic_solutions_hk_0/solution_0000074.mat")

nz = mesh_params.nz; 
mesh_params.w0 = u(1:nz); 
mesh_params.w0z = mesh_params.Dz*u(1:nz); 
w0 = u(1:nz); 
w0k = u(nz+1:2*nz);  
k = p(5); 
b = p(3); 
c = 0; 

u0 = [w0; w0k; 15*w0k; c; c; c; b; k]; %factor of gives good initial guess for fast convergence. 

%converge initial guess 
my_rhs2 = @(u) SHE357_1D_dHdk_0_dHdkk_0(u,p,mesh_params,0);
options = optimset('Jacobian','off','Display','iter','MaxIter',5e10,'MaxFunEval',5e10,...
    'Algorithm','levenberg-marquardt');
% options = optimset('Jacobian','off','Display','iter','MaxIter',500);
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs2,u0,options);

%% Continue in the bifurcation parameter mu using Daniele's continuation code
% Various function handles
problemHandle            = @(u,p)SHE357_1D_dHdk_0_dHdkk_0(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_SHE357_dHdk_0_dHdkk_0(step,u,p,mesh_params);
computeEigenvaluesHandle = @(u,p) ComputeSpectrum_hkk_0(u,p,mesh_params);
plotSpetcrumHandle       = @(d,p,parentHandle) PlotSpectrum(d,p,parentHandle); 
stepperPars.iContPar      = 2; % parameter to continue 
stepperPars.s0            = -1; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = 10; 
stepperPars.pMin          = -50.0;
stepperPars.pMax          = 50;
stepperPars.maxSteps      = 1e3;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-3;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',30,'TolFun',1e-12); 
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = 'codim2_pt';
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpetcrumHandle;  
 % stepperPars.PlotSpectrum = []; 
 stepperPars.PlotSolution  = []; 
stepperPars.PlotBranchVariableId = 7;  
PLbranch = SecantContinuation(problemHandle,u_out,p,stepperPars);
