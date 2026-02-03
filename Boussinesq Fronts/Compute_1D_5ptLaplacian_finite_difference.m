%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 1D finite-difference Laplacian - Neumann bcs on r=[0,Lrh]
% Inputs: Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
% Note: returned matrices are sparse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,L,Dx,w,D4x] = Compute_1D_5ptLaplacian_finite_difference(Nrh,Lrh)

hx = Lrh/(Nrh-1);

r = (0:Nrh-1)'*hx;r(1) = 1;                      % radial mesh
R = sparse(1:Nrh,[1:Nrh],1./r,Nrh,Nrh); r(1) = 0;

%% Finite-difference matrices
ex = ones(Nrh,1);

Ix = sparse(1:Nrh,[1:Nrh],ex,Nrh,Nrh); % I

Dx = spdiags([ex -8*ex 0*ex 8*ex -ex],-2:2, Nrh, Nrh);
Dx(1,:)   = 0; Dx(2,2)   = 1;
Dx(Nrh,:) = 0; Dx(Nrh-1,Nrh-1) = -1;
Dx = Dx/(12*hx);

D2x = spdiags([-ex 16*ex -30.*ex 16*ex -ex], -2:2, Nrh, Nrh);
D2x(1,2)= 32; D2x(1,3)= -2;
D2x(2,1)= 16; D2x(2,2)= -31; 
D2x(Nrh,:) = D2x(1,Nrh:-1:1);
D2x(Nrh-1,:) = D2x(2,Nrh:-1:1);

D2x = D2x/(12*hx^2);

D4x = spdiags([-ex 12*ex -39*ex 56*ex -39*ex +12*ex -ex],-3:3, Nrh, Nrh);
D4x(1,2) = -78; D4x(1,3) = 24; D4x(1,4) = -2;
D4x(2,2) =  68; D4x(2,3) = -40;
D4x(3,2) = -40;
D4x(Nrh,:) = D4x(1,Nrh:-1:1);
D4x(Nrh-1,:) = D4x(2,Nrh:-1:1);
D4x(Nrh-2,:) = D4x(3,Nrh:-1:1);

D4x = D4x/(6*hx^4);


%% 1D
L = D2x;

%% 2D radial
% L = D2x + R*Dx; % d_rr + d_r/r
% L(1,1) = 2*D2x(1,1); L(1,2) = 2*D2x(1,2); % 2d_rr at r= 0 

%% 3D radial
% L = D2x + 2*R*Dx; % d_rr + 2*d_r/r
% L(1,1) = 3*D2x(1,1) + 1; L(1,2) = 3*D2x(1,2); % 3d_rr at r= 0
if (nargout>2)
  % w = [1, 2*ones(1,Nrh-2), 1]*hx;   % trapezoidal weights for intergration int = w*u
    w = [1, 2*ones(1,Nrh-2)+2*mod([1:Nrh-2],2),1]*hx/3; % simpson's rule weights for integration int = w*u
end
