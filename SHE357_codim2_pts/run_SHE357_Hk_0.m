%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes 2-pi periodic stationary solutions to the stationary 
% cubic quintic septic Swift-Hohenberg equation (SH357)
%         −(1 + k^2*d_zz)^2[u] + μ*u - a u^3 + bu^5 - u^7  = 0 
% with the additional constraint H_k = 0, where H is the spatial
% Hamiltonian. Initial guess is obtained from the output from the script
% run_SHE357_periodic.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load solution from previous computations and the computational mesh 
load('periodic_solutions/mesh_params.mat'); 
load("periodic_solutions/solution_0000033.mat")

mesh_params.w0 = u(1:end-1); 
mesh_params.w0z = mesh_params.Dz*u(1:end-1); 
w0 = u(1:end-1); 
k = p(5); 
c = 0; 

%output vector u0 = [u,uk,phase1,phase2,k]
u0 = [w0; w0; c; c; k];

% converge an initial guess--modify as needed 
disp('Converging periodic orbit up');
my_rhs = @(u) SHE357_1D_dHdk_0(u,p,mesh_params,0);
options = optimset('Jacobian','off','Display','iter','MaxIter',5e3);
% options = optimset('Jacobian','off','Display','iter','MaxIter',500);
[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,u0,options);

%% Continue in the bifurcation parameter mu using Daniele's continuation code
% Various function handles
problemHandle            = @(u,p)SHE357_1D_dHdk_0(u,p,mesh_params);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_cubic(u,p,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_SHE357_dH_dk_0(step,u,p,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpetcrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle); 
stepperPars.iContPar      = 3; % parameter to continue 
stepperPars.s0            = 0.01; % minus sign here steps backwards in parameter
stepperPars.sMin          = 1e-2;
stepperPars.sMax          = 0.1;
stepperPars.pMin          = -50.0;
stepperPars.pMax          = 100;
stepperPars.maxSteps      = 250;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',5000);
stepperPars.optNonlinIter = 10;

stepperPars.dataFolder    =  'periodic_solutions_hk_0'; 
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = [];   
stepperPars.PlotBranchVariableId = 8; 
branch = SecantContinuation(problemHandle,u_out,p,stepperPars);

%save mesh parameters. If the continuations is stopped prematurely, this
%may have to be done manually. 
if ~exist(stepperPars.dataFolder, 'dir')
    mkdir(stepperPars.dataFolder)
end
save(fullfile(stepperPars.dataFolder, 'mesh_params'), 'mesh_params');
pause(0.01)

%generate figure 11: 
subplot(131); plot(branch(:,3),branch(:,end)); 
subplot(132); plot(branch(:,3),branch(:,end-1)); 
subplot(133); plot(branch(:,3),branch(:,8)); 