%% Elfouhaily Andreas Divinyi 2018 Chalmers 
% Implemented from equation found at http://www.oceanopticsbook.info/view/surfaces/level_2//wave_variance_spectra__examples
% Input windspeed,development,degreees and wavenumber
% Degrees work in progress for spread function
function s=elfouhaily(U_10,omegac,phi,k)
%Constants
alpha = 0.0081; % Enhet? [] och betydelse
beta= 1.25;     % Enhet? []
g=9.82;         %[m/s^2] gravitional acceleration
%U_10=5:5:15;          %[m/s] windspeed 19.5m above seasurface
%omegac=1;       % Depends on how developed the sea is
Cd10N=0.00144;  % Drag
ustar=sqrt(Cd10N).*U_10(1);
a0=0.1733;
ap=4;
km=370;         %[rad/m]
cm=0.23;        %[m/s]
am=0.13*ustar./cm;
%Beror på omegac större än 1 lägg till 6*log10(omegac) This part will have
%to be moved to 
if omegac <=1
    gamma=1.7;
else
    gamma=1.7+6*log10(omegac);
end

sigma=0.08*(1+4*(omegac^-3));
alphap=0.006*omegac^(0.55);

% There's a if-else dependency on cm for ustar
if ustar <= cm
    alpham=0.01*(1+log(ustar./cm));
else
    alpham=0.01*(1+3*log(ustar./cm));
end

k0=g./(U_10(1).^2);
kp=k0*(omegac^2);
cp=sqrt(g./kp);
c=sqrt((g./k).*(1+(k/km).^2));

% Creates spectrum
Lpm=@(k)exp(-1.25.*((kp./k).^2));
Gamma=@(k)exp((-1./(2.*(sigma.^2))).*((sqrt(k/kp)-1).^2));
Jp=@(k)gamma.^Gamma(k);
Fp=@(k)Lpm(k).*Jp(k).*exp(-0.3162*omegac.*(sqrt(k/kp)-1));
Fm=@(k)Lpm(k).*Jp(k).*exp(-0.25*(((k./km)-1).^2));
Bl=@(k)0.5*alphap.*(cp./c).*Fp(k);
Bh=@(k)0.5.*alpham.*(cm./c).*Fm(k);

S=@(k)(Bl(k)+Bh(k))./(k.^3);
% Phi=@(k,phi)1/(2*pi)*(1+tanh(a0+ap.*(c./cp).^(2.5)+am*(cm./c).^(2.5)).*cos(phi*2));

% Psi=@(k,phi)S(k).*Phi(k,phi)./k;


s=S(k);