function [x,y] = compute_leaf(uInput,p,phi,tf,epsilon,phase,manifLength)

%inputs
%phi : Mondoromy matrix
%t0 :

c = -1; %for cubic SH
k = p(2); %wavenumber
sigma = p(4);
g = p(5);
nn = p(6);

N = length(uInput);
q = [-N/2:N/2-1]'*k;
uHat = fftshift(fft(uInput(:,1))/N);

%interpolate the solution & derivatives using the Fourier series
uHat = uHat;
u = @(z) (real(sum(uHat.*exp(1i*q*(z)))));
u2 = @(z) (real(sum(1i*q.*uHat.*exp(1i*q*(z)))));
u3 = @(z) (real(sum((1i*q).^2.*uHat.*exp(1i*q*(z)))));
u4 = @(z) (real(sum((1i*q).^3.*uHat.*exp(1i*q*(z)))));

[V,D] = eig(reshape(deval(phi,2*pi/k),[4,4]),'vector');

if tf > 0 %unstable manifold (integrating forward in time)
    [~,ind] = max(abs(D));
    Y = real(V(:,ind)); %get eigenvector of monodromy matrix
    Y = Y/norm(Y)*sign(Y(3));
else %perturbation on the stable manifold
    [~,ind] = min(abs(D));
    Y =real( V(:,ind)); %get eigenvector of monodromy matrix
    Y = Y/norm(Y)*sign(Y(3));
end

%build vector of perturbations

pert = reshape(deval(phi,phase),[4,4])*Y;
pert = pert/norm(pert,2);
X0= [u(phase);u2(phase);u3(phase);u4(phase)] + epsilon*pert;

opts = odeset('RelTol',100*eps,'AbsTol',eps...
    ,'NormControl','off','InitialStep',eps...
    ,'Events',@(t,y)exitEvent(t,y),'Refine', 1);
[x,y] = ode113(@(t,y)fun(t,y,c,sigma,g,nn),[linspace(0,tf,manifLength)],[X0],opts);


    function dydx = fun(x,y,c,sigma,g,nn)

        dydx = [y(2,:)
            y(3,:)
            y(4,:)
            c*y(1,:) - nn*y(1,:).^2 - g*y(1,:).^3  - sigma*y(3,:)];
        dydx = dydx(:);
    end

    function [position, isterminal, direction] = exitEvent(t, y)
        %if the ODE solution gets too large, halt integration
        position = y(1); 
        isterminal = 1;  % halt integration
        direction = 0;
    end
end





