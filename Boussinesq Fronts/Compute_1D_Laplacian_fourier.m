%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 1D fourier Laplacian - r=[-Lrh,Lrh]
% Inputs: 2*Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,L,Dx,w] = Compute_1D_Laplacian_fourier(Nrh,Lrh)

% Fourier differentiation matrix first order for y between -pi and pi
h = 2*pi/Nrh;  t = h*(1:Nrh); r = Lrh*(t-pi)/pi;
column = [0 .5*(-1).^(1:Nrh-1).*cot((1:Nrh-1)*h/2)]';
Dt = toeplitz(column,column([1 Nrh:-1:2]));
Dtt = toeplitz([-pi^2/(3*h^2)-1/6 .5*(-1).^(2:Nrh)./sin(h*(1:Nrh-1)/2).^2]);


%% 1D
Dx = (pi/Lrh)*Dt;
L  = (pi/Lrh)^2*Dtt;

if (nargout>2)
   w = 2*Lrh*ones(1,Nrh)/Nrh;   % trapezoidal weights for intergration int = w*u
end
