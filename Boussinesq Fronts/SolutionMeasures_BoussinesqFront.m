function F = SolutionMeasures_BoussinesqFront(step,u,p,mesh_params)

global ur0_left ur0_right;

% Solution measures (they are displayed on screen)
% by setting stepperPars.PlotBranchVariableId = k,
% the kth solution measure is plotted in the
% bifurcation diagram on the fly

A  = p(1);
R  = p(2);
B  = p(3);
D1 = p(4);
D2 = (A/R)^2;
m  = p(5);
d  = p(6);
phi= p(7);

kxm = u(end-1); 
kxp = u(end);

F = [kxm,kxp]; %output far-field wavenumbers
end
