function [phi] = get_monodromy_cubic(uInput,p)

c = p(1); %speed
k = p(2); %wavenumber
A = p(3); %constant of integration (usually 0);
alpha = p(4);
g = p(5);
nn = p(6); 

t0 = 0;
tf = 2*pi/k;
uSol = uInput;

dk = k; N = length(uInput);
q = [-N/2:N/2-1]'*dk;
uHat = fftshift(fft(uInput(:,1))/N);
q = [q;-q(1)];

uHat = [uHat;uHat(1)];

%interpolant based on Fourier series
u = @(z) (real(sum(uHat.*exp(1i*q*(z)))));
u2 = @(z) (real(sum(1i*q.*uHat.*exp(1i*q*(z)))));
u3 = @(z) (real(sum((1i*q).^2.*uHat.*exp(1i*q*(z)))));
u4 = @(z) (real(sum((1i*q).^3.*uHat.*exp(1i*q*(z)))));

%compute the transition matrix for kdv5 equation
phi0 = reshape(eye(4),[],1);

opts = odeset('RelTol',100*eps,'AbsTol',eps...
    ,'NormControl','off','InitialStep',eps);


phi = ode113(@(z,phi)dphidx(z,phi,@(z)u(z),c,alpha,A,g,nn),[t0,tf],phi0,opts);
 
if min(abs(eig(reshape(deval(phi,2*pi/p(2)),[4,4]))'-1)) > 1e-3
    error('Monodromy matrix does not have repeated eigenvalues at +1. This is likely due to an inaccurate periodic solution')
end


    function phi_prime = dphidx(t,phi,u,c,alpha,A,g,nn)
        phi = reshape(phi,4,4);
        A = [0                                    1   0      0;
             0                                    0   1      0;
             0                                    0   0      1;
            -3*g*u(t).^2 + c - 2*nn*u(t)+ 0*A     0  -alpha  0];
        
        phi_prime = reshape(A*phi,[],1);

    end


end