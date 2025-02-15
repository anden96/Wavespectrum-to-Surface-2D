%% spec2surf

clc
clear all

g = 9.81;
L = 100;
N = 64;
u = linspace(-((N/2)-1),N/2,L);
uv = u.*u';
w = linspace(-0.1,50*pi,L);
k = (w.^2)/g;
a = zeros(L,L);

Psi1s = elfunW(k,1);

%Check ishermitian(z) to see that it returns 1
    rho1 = normrnd(0,1,L,L);
    sigma1 = normrnd(0,1,L,L);
    Psi2s = circshift(Psi1s,-51);
    zhat01 = (1./sqrt(2)).*(rho1 + j.*sigma1).*sqrt(Psi2s); 

    z = (1./sqrt(2)).*(zhat01+ctranspose(zhat01));
    
    Z = ifft2(z);
    imag(Z)
    Y = real(Z);

ishermitian(z)

surf(Y)
