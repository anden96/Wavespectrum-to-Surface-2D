 function E2D=ECKV2D_k_phi(kpos,phirad,U10)

% ; Given positive angular wavenumbers kpos in rad/m and angle phi in radians, and wind speed U10 in m/s, 
% ; this routine returns the 2D directional gravity-capillary energy spectrum Psi(k,phi) of 
% ; Elfouhaily, T., B. Chapron, K. Katsaros, and D. Vandmark, 1997. A unified directional spectrum for 
% ; long and short wind-driven waves. J Geophys Res, vol 102(C7), 15781-15796.
% ; Psi(k,phi) is given by ECKV Eq. (67) and related equations.
% ; This is a one-sided energy spectrum valid for positive wavenumbers in the full gravity-capillary range.
% 
% ; Input:
% ;   k = the angular wavenumber k, in rad/m
% ;   phi = the angle in rad, -pi to pi or 0 to 2pi 
% ;        k and phi can be single wavenumbers or an array of positive wavenumbers
% ;   U10 = the wind speed at 10 m above sea level
% ;         Note: U10 must be passed via a common block if the function is called by a numerical
% ;            integration routine such as qromb, which requires as input a function routine 
% ;            with just one argument for the abscissa
% 
% ; Output:
% ;   Psi(k,phi) = the one-sided energy spectrum for the given (k,phi) values and U10, 
% ;                in units of m^2/(rad/m)^2 (shown as m^4/rad^2 in ECKV Table A1)
% ;                This E(kx,ky) is the Psi(kx,ky) = k Psi(k,phi) where Psi(k,phi) is given by Eq. (67)
%  
% ;------ parameters defining the sea state

g = 9.82;

% ; Omegac = 0.84 for a "fully developed" sea (corresponds to Pierson-Moskowicz)
% ;        = 1 for a "mature" sea (used in ECKV Fig 8a)
% ;        = 2 to 5 for a "young" sea; max allowed value is 5
Omegac = 0.84

% ;  convert U10 = wind at 10 m to the friction velocity u* using ECKV Eq. 61
Cd10N = 0.00144; %          ; value deduced from ECKV Fig 11
ustar = sqrt(Cd10N)*U10;

% ; values from ECKV Eq. 59
ao = 0.1733;
ap = 4.0;
cm = 0.23;
am = 0.13*ustar/cm;
km = 370.0; %     ; corresponds to lambda = 1.7 cm for gravity-capillary wave boundary

if Omegac <= 1  
   gamma = 1.7;
else
   gamma = 1.7 + 6.0*log10(Omegac);
end

sigma = 0.08*(1.0 + 4.0*Omegac^(-3));
twosig2 = 1.0/(2.0*sigma*sigma);
alphap = 0.006*Omegac^(0.55); % ECKV Eq. 34

% ; ECKV Eq. 44:
if ustar <= cm 
   alpham = 0.01*(1.0 + log(ustar/cm));
else
   alpham = 0.01*(1.0 + 3.0*alog(ustar/cm));
end

% ;  get the wavelength corresponding to the max of the wavenumber spectrum
ko = g/(U10*U10);
kp = ko*Omegac*Omegac;
wavemax = 2.0*pi/kp;

cp = sqrt(g/kp);
Capomega = Omegac; %  ; = U10/cp for thetabar = 0

% ;------ now evaluate the energy spectrum equation

Nk = length(kpos);
% ;phirad = phideg/!radeg
Nphi = length(phirad);
E2D = zeros(Nk,Nphi);

for ku=1:Nk 
  for phiv=1:Nphi
    k = kpos(ku);
    phi = phirad(phiv);
    
%     ;  define the S(k) spectrum for the given k value
    
    c = sqrt((g/k)*(1.0 + (k/km)^2)); %; = sqrt(g/k + Tw*k/rhow)
    sqrtkkp = sqrt(k/kp);
    
%     ;  the low-frequency (long-wave) part  B_l of Eq. (31)
    Capgamma = exp(-twosig2*(sqrtkkp - 1.0)^2); %; below Eq. (3)
    fJp = gamma^Capgamma; % Eq. (3)
    fLpm = exp(-1.25*(kp/k)^2) %; Eq. (2)
    Fp = fLpm*fJp*exp(-0.3162*Capomega*(sqrtkkp - 1.0)) %; Eq. (32)
    Bl = 0.5*alphap*(cp/c)*Fp;    %; sqrt(k/kp) = cp/c
%     ;Slow = Bl/k^3
    
%     ; the high-frequency (short-wave) part B_h of Eq. (40)
%     ; NB: A typo in Eq. (41) of the Elfouhaily paper omits the Lpm*Jp factor
%     ; from the definition of Fm; compare with Fp in Eq. (32)
    Fm = fLpm*fJp*exp(-0.25*(k/km - 1.0)^2); % Eq. (41) with typo fixed
    Bh = 0.5*alpham*(cm/c)*Fm; % Eq. (40)
%     ;Shigh = Bh/k^3
    
    Deltak = tanh( ao + ap*(c/cp)^(2.5) + am*(cm/c)^(2.5) ) % ; Eq. (57)
    
%     ; combine as in Eq. (67) to get (Note: this is Psi(k,phi), so 1/k^4)
    
%     ; normal runs return the energy spectrum Psi(k,phi)
%     ;   E2D[ku,phiv] = (0.5/!pi)*(Bl + Bh)/k^4 * (1.0 + Deltak*cos(2.0*phi))
%     ; return only downwind values and double the spreading function
    if phi >= -0.5*pi && phi <= 0.5*pi
      E2D(ku,phiv) = (1.0/pi)*(Bl + Bh)/k^4 * (1.0 + Deltak*cos(2.0*phi));
    else
      E2D(ku,phiv) = 0.0;
    end  
    
%     ; use this to return just the core of the spreading function seen in Eq. (49)
%     ; use this for reproducting Fig. 10a
%     ;E2D[ku,phiv] = (1.0 + Deltak*cos(2.0*phi))

    end %; phiv
  end %; ku
<<<<<<< HEAD
  
  

=======
>>>>>>> refs/remotes/origin/Beta-branch
