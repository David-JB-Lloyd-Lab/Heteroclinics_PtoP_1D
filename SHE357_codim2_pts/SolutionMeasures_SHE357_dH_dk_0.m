function F = SolutionMeasures_SHE357_dH_dk_0(~,uu,p,mesh_params)

mu = p(1);
a = p(2);
b  = p(3);
sigma = p(4);
% k  = p(5);
gamma = p(6); 

% Auxiliary variables
n = mesh_params.nz;
u = uu(1:n);
uk= uu(1+n:2*n);
c1 = uu(2*n+1);
c2 = uu(2*n+2);
k  = uu(2*n+3);

w = mesh_params.wz;
D2z = mesh_params.D2z;
Dz = mesh_params.Dz;

I = speye(n);
% E = w*(0.5*((I + k^2*D2z)*u).^2 - 0.5*mu*u.^2 - nu*u.^3/3 + g*u.^4/4);
E = w*((D2z*u).^2 - sigma*(Dz*u).^2/2 - 0.5*mu*u.^2 + a*u.^4/4 - b*u.^6/6 + u.^8/8);
H = -0.5*k^4*(D2z*u).^2 + sigma*k^2*(Dz*u).^2/2 + (1-mu)*u.^2/2  ...
    + 1/4*a*u.^4 - b*u.^6/6 + u.^8/8 + k^4*(Dz*u).*(Dz*D2z*u);

A = mean(sigma*k*(Dz*u).^2 - 2*k^3*(D2z*u).^2);

Ak= mean(sigma*(Dz*u).^2 + 2*sigma*k.*(Dz*u).*(Dz*uk) ...
    - 6*k^2*(D2z*u).^2 - 4*k.^3*(D2z*u).*(D2z*uk));


%adding some bits to compute the parametric derivatives of u
opts = optimset('Display','off','Algorithm','Levenberg-Marquardt',...
                                     'DerivativeCheck','off',...
                                     'Jacobian','off',...
                                     'MaxIter',15);

% uk = fsolve(@(uk)SHE357_1D_get_Hkk(uk,[u;k],p,mesh_params),u,opts);
% 
dHdk = -2*k^3*(D2z*u).^2 -k^4*(D2z*u).*(D2z*uk) ...
    + sigma*k*(Dz*u).^2 + sigma*k^2*(Dz*u).*(Dz*uk) ...
    + (1-mu)*u.*uk + a*u.^3.*uk - b*u.^5.*uk + u.^7.*uk ...
    + 4*k^3*(Dz*u).*(Dz*D2z*u) + k^4*(Dz*uk).*(Dz*D2z*u) ...
    + k^4*(Dz*u).*(Dz*D2z*uk);


F = [mesh_params.wz*(u.^2)/mesh_params.Lz, ...
    E,...
    H(1),...
    A,...
    dHdk(1),...
    k,...
    Ak, ... 
    norm(u,2),...
    k]; 
end
