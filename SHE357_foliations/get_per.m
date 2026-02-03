function [U,H] = get_per(u,p,opts)

c = p(1); %speed
k = p(2); %wavenumber
alpha = p(4);
g = p(5);
nn = p(6);

N = length(u);
q = fftshift(-N/2:N/2-1)';
z = pi/N*[-N/2:N/2-1]; 

F = @(u) -c*u + nn*u.^2 + g*u.^3 + ifft((-alpha*(q*k).^2 + (q*k).^4).*real(fft(u)));
U = fsolve(F,u,opts);

if max(U)-min(U) < 1e-8 %check if the zero solution was computed.
    warning('solution went to a constant -- modify initial guess..')
end

%definitions to compute the Hamiltonian. 
Up = ifft(1i*q.*fft(U),'symmetric');
Upp = ifft(1i*q.*fft(Up),'symmetric');
Uppp = ifft(1i*q.*fft(Upp),'symmetric');
U4p = ifft(1i*q.*fft(Uppp),'symmetric');

H = (-c/2*U.^2 + nn*(1/3)*U.^3 + (g/4)*U.^4 + 0.5*alpha*k.^2*(Up).^2 + k.^4*Uppp.*Up ...
    - k.^4*0.5*(Upp).^2); %hamiltonian value.
end