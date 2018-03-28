% @ECKV2D_k_phi ; the one-sided, 2D ECKV spectrum with the symmetric ECKV Phi(k,phi) spreading function
% NEED TO IMPLEMENT ECKV2D_k_phi in matlab to 
% pro gen2Dsurf

% ; This routine shows how to use a one-sided 2D variance spectrum Psi_1s(k,phi) and 2D inverse FFTs 
% ; to generate a random 2D sea surface.
% ; (To generate a time series of 2D surfaces for creation of videos, use routine cgAnimate_2D_SeaSurface.pro)
% 
% ; This routine created Tutorial Fig. 3.3
% 
% ; The Elfouhaily et al (1997) directional variance spectrum is used (the Psi(k,phi) of their Eq. 67). 
% ; This is computed by routine ECKV2D_k_phi.
% 
% ; This version uses angular frequencies kx,ky = 2 pi/wavelength [rad/m] as the spatial frequency variables.
% 
% ; Special care is taken to define the values at the ky and ky Nyquist frequencies, which are
% ; always special cases.
% 
% ; Results are shown as 4 contour plots:
% ; The upper left  panel is the one-sided wave variance spectrum Pis_1S(kx,ky)
% ; The upper right panel is the real part of the random zhat(kx,ky)
% ; The lower left  panel is the imaginary part of the random zhat(kx,ky)
% ; The lower right panel is the generated random sea surface z(x,y)
% 
% ; ***************************************************************************************************
% ; * This code is copyright (c) 2016 by Curtis D. Mobley.                                            *
% ; * Permission is hereby given to reproduce and use this code for non-commercial academic research, *
% ; * provided that the user suitably acknowledges Curtis D. Mobley in any presentations, reports,    *
% ; * publications, or other works that make use of the code or its output.  Depending on the extent  *
% ; * of use of the code or its outputs, suitable acknowledgement can range from a footnote to offer  *
% ; * of coauthorship.  Further questions can be directed to curtis.mobley@sequoiasci.com.            *
% ; ***************************************************************************************************
%   * Translated from IDL to MATLAB code by Andreas Divinyi 2018 during a
%   bachelor project at Chalmers, Gothenburg, Sweden
% ;-----------------------------------------------------------------------------------------------------------
% ;| This routine uses Coyote Graphics routines rather than standard IDL for plotting; see www.idlcoyote.com |
% ;| and "Coyote's Guide to Traditional IDL Graphics: Using Familiar Tools Creatively" by David Fanning      |
% ;| Any comments below like "CG page xxx" refer to this book.                                               |
% ;| CG routines are wrappers for standard IDL routines, which make it easier to do certain things like      |
% ;| create vector eps output files; cgPlot replaces the standard Plot, etc.                                 |
% ;-----------------------------------------------------------------------------------------------------------
%  Use of MATLAB standard graphics packets
% 
% ; Written 13 August 2014 by Curtis Mobley and modified by Andreas Divinyi
% 2018
% 
% 
% ;-----------------  BEGIN USER INPUT  ----------------------

idbug=0; % Set to 1 for debugging output

% ; ***** DEFINE THE PHYSICAL REGION AND SAMPLING ***** 
% ; define the wind speed at 10 m above MSL for use in the variance spectrum
U10 = 5.0; % [m/s]

Nx = 64; % number of samples of sea surface elevation to be generated in the x direction; MUST be a power of 2 for the FFT
Ny = 64; % number of samples of sea surface elevation to be generated in the y direction; MUST be a power of 2 for the FFT

InputOption=2;

switch InputOption
    case 1
        labelInOpt = ' Spatial resolution is defined by the sample spacing and number of sample points';
        % smallest resolvable wavelength is 2*Deltax
        Deltax = 0.1; %0.0025
        Deltay = 0.1;
        Lx = Deltax*Nx;
        Ly = Deltax*Ny;
    case 2
          labelInOpt = ' Spatial resolution is defined by the spatial region and number of sample points';
        % largest resolvable wave is the length of the spatial region
        Lx = 100.0; %length of surface region in the x direction, in meters 
        Ly = 100.0; %length of surface region in the y direction, in meters 
        Deltax = Lx/Nx;
        Deltay = Ly/Ny;
end

% ; Define root names for output files.
% ; File names will have the form root name + wind speed + x spatial size + x sampling, e.g.
% ; PlotFileName = 'Fig3.3_U10_L100_N1024.eps'
% ; The info on wind speed, spatial size, and sampling size will be appended below.

PlotRootName = 'Fig3.3_ECKV';

% ; Optionally save a text file z2D_PlotRootName.txt of (x,y,z) values for postprocessing
% ; These z(x,y) files can be plotted by routines cgPlot2Dsurf_3D and cgPlot2Dsurf_contour
% ; Note: these files can be very large for large Nx, Ny
isave = 1; % = 1 to save output, 0 to not save

% Seed probably not needed for matlab
% ; define an initial seed for random number generation
% seed =  833202L; 538832L ;33202L; 838392L ;33202L;
% ; Note: seed = 833202L was used on a 32 bit computer to generate Tutorial Fig. 3.3.
% ; Use of a different seed, or this seed on a 64 bit computer, will generate a different
% ; surface realization because a different sequence of random numbers will be generated.

% ;-----------------  END USER INPUT  ----------------------
% 
% ; create output file names
% Create a Ustring 
Ustring=num2str(U10);
if U10 < 10
    Ustring=Ustring(1);
elseif U10 >= 10
    Ustring=[Ustring(1) Ustring(2)];
end

Lstring=num2str(Lx);
if Lx < 10
    Lstring=Lstring(1);
elseif Lx >=10 && Lx < 100
    Lstring = [Lstring(1) Lstring(2)];
elseif Lx >= 100
    Lstring= [Lstring(1) Lstring(2) Lstring(3)];
end

Nstring = num2str(Nx);
if Nx < 10
    Nstring = Nstring(1);
elseif Nx >=10 && Nx < 100
    Nstring=[Nstring(1) Nstring(2)];
elseif Nx >= 100 && Nx <1000
    Nstring=[Nstring(1) Nstring(2) Nstring(3)];
elseif Nx >= 1000
    Nstring=[Nstring(1) Nstring(2) Nstring(3) Nstring(4)];
end

PlotFileName=[PlotRootName '_U' Ustring '_L' Lstring '_N' Nstring '.eps'];

% ; ***** COMPUTE THE SPATIAL FREQUENCIES *****
kxmin = 2.0*pi/Lx; % min frequency (min spatial frequency); frequency of longest resolvable wave
kymin = 2.0*pi/Ly;
minxwave = 2.0*Deltax; % min resolvable wavelength (2 point wave)
minywave = 2.0*Deltay; % min resolvable wavelength (2 point wave)
Nyquistx = pi/Deltax;  % Nyquist spatial freq in rad/m
Nyquisty = pi/Deltay;  %Nyquist spatial freq in rad/m
kxmax = Nyquistx;
kymax = Nyquisty;

% Print statements -NEEDS fixing from print to disp
%{
print,' '
print,' Generation of a 2D sea surface from a one-sided, 2D variance spectrum Psi1s(kx,ky) = Psi(k,psi)'
print,' '
print,' Note: This routine does NOT account for unresolved variance'
print,' '
print, format='(4x,a,2i5)',   'Numbers of spatial samples = Nx, Ny =',Nx,Ny
print, format='(4x,a,2f7.4)', 'Sample spacings Deltax, Deltay [m] =', Deltax,Deltay
print, format='(4x,a,2f8.1)', 'Max resolvable wavelengths = Lx, Ly [m] =',Lx,Ly
print, format='(4x,a,2f7.4)', 'Min resolvable wavelengths = 2*Deltax, 2*Deltay [m] =', minxwave,minywave
print, format='(4x,a,2f7.2)', 'Min resolvable frequencies kxmin, kymin [rad/m] = ',kxmin ,kymin
print, format='(4x,a,2f7.2)', 'Max resolvable (Nyquist) frequencies kxmax, kymax [rad/m] = ',kxmax ,kymax
print, format='(4x,a,i12)',   'Seed for random number generation =',seed
print,' '
%}

% ; ***** NOTE ON FREQUENCY ARRAY ORDER:
% 
% ; When plotting 1-sided power spectra, it is customary to plot only positive frequencies in the order
% ; kpos[1 to N/2] = [1, 2,..., N/2-1, N/2(Nyquist freq)]*(2pi/L)  ; total of N/2 pos frequencies
% 
% ; When plotting 2-sided power spectra, and when going from z(x,y) to zhat(kx,ky) = FFT2D{z(x,y)},
% ; mathematicians usually write the frequencies in the order
% ; kmath[1 to N] = [-(N/2-1),...,-2,-1,0,+1,+2,...,+N/2-1,+N/2(Nyquist freq)]*(2pi/L)  ; total of N freqs in rad/m
%  
% ; IDL (and Matlab) stores these frequencies in the order (see the IDL documentation)
% ; kIDL[1 to N] = [0.0, 1 2, ..., (N-1)/2, N/2(Nyquist freq), -(N-1)/2, ..., -2, -1]*(2pi/L) ; total of N freqs
% ; This is the (kx,ky) order of zhat(kx,ky) arrays returned by IDL routine FFT.
% 
% ; The FFT order is generated as follows:
% ; FFT = FINDGEN((Nx-1)/2) + 1.0 ; = [1,2,...,(N-1)/2]
% ; kFFT = (2.0*!pi/L)*[0.0, FFT, N/2.0 ,-N/2.0+FFT] ; in rad/m
% 
% ; The math order can be obtained from the FFT order by a circular shift by N/2-1 to the right:
% ; kmath = SHIFT(kFFT,Nx/2-1)
% 
% ; The FFT order can be obtained from the math order by a circular shift by N/2-1 to the left:
% ; kFFT = SHIFT(kmath,-(N/2-1))
% 
% ; When defining the one-sided variance spectrum Psi1s for random sampling, we will use (kvec,phi) with 
% ; vector kvec = sqrt(kx^2 + ky^2) including zero frequencies and -pi/2 <= phi <= pi/2.  This 
% ; corresponds to 0 and positive kx, and positive and negative ky.  These frequencies are 
% ; denoted by kx1S = [0.0, kxpos] (Nx/2 +1 values) and ky1S = kymath (Ny values)
% 
% ; These forms of frequency arrays are computed next for use below.
% 
% ;----- Compute positive frequencies used for plotting the 1-sided variance spectrum; 
% ;   these do not include 0 frequency;
% ;   even spacing in frequency k; 
% ;   the last frequency is the Nyquist frequency
Nkxpos = Nx/2;
Nkypos = Ny/2;
% ; spacing for positive frequencies; Delta includes the 2pi/L values via the kmax and kmin values
Deltakx = (kxmax - kxmin)/(Nkxpos-1); %-1 relevant in matlab?  
Deltaky = (kymax - kymin)/(Nkypos-1);
kxpos = kxmin + Deltakx*(1:Nkxpos); % in rad/m
kypos = kymin + Deltaky*(1:Nkypos);

% ;----- Compute the FFT frequencies needed by the FFT routine:
xFFT = 1:((Nx-1)/2); % = [1,2,...,(Nx-1)/2]
yFFT = 1:((Ny-1)/2);
kxFFT = (2.0*pi/Lx)*[0.0, xFFT, Nx/2.0 ,-Nx/2.0+xFFT]; % in rad/m
kyFFT = (2.0*pi/Ly)*[0.0, yFFT ,Ny/2.0, -Ny/2.0+yFFT];
%To plot the kyFFT and kxFFT - Looks OK 
% figure(1)
% subplot(1,2,1)
% plot(kyFFT)
% subplot(1,2,2)
% plot(kxFFT) 

% ;----- Compute the math frequencies
% Worth controlling if any differences with SHIFT and circshift
% Confirmed that circshift does the same as shift, ref: Harris Geospatial
% Solutions, IDL reference. 18-03-04
kxmath = circshift(kxFFT,Nx/2-1);
kymath = circshift(kyFFT,Ny/2-1);
%To plot the kymath and kxmath - Looks OK  
% figure(2)
%  subplot(1,2,1)
%  plot(kymath)
%  subplot(1,2,2)
%  plot(kxmath) 

% ;----- Compute the 1S frequencies

kx1S = [0, kxpos] ;
ky1S = kymath;

%To plot the kx1S and ky1S  
figure(3)
title('kx1S and ky1s')
subplot(1,2,1)
plot(kx1S)
subplot(1,2,2)
plot(ky1S) 


% ;if(idebug eq 1) then begin
% ;  print,' kxpos  = ',kxpos
% ;  print,' kxFFT  = ',kxFFT
% ;  print,' kxmath = ',kxmath
% ;  print,' kx1S   = ',kx1S
% ;  print,' kypos  = ',kypos
% ;  print,' kyFFT  = ',kyFFT
% ;  print,' kymath = ',kymath 
% ;  print,' ky1S   = ',ky1S 
% ;endif
 

% ; ***** SET UP THE PLOT ***** 
% MOST OF THIS WILL PROBABLY BE HANDLED AUTOMATICALLY BY MATLAB
% ; Set up a multipanel layout
% ;   (Note: do NOT use /noerase with multipanel plots!) - Along the lines of hold on in matlab


nplotcol = 2; % number of plot columns
nplotrow = 2; % number of plot rows

% !p.multi = [0,nplotcol,nplotrow] 
% 
% ; Set up for a vector eps output file; a high-quality png file will also be created
% 
% !x.margin = [5,7]
% !y.margin = [4,3]
% !x.Omargin = [2,2]
% !y.Omargin = [1,3]
% 
% aspectratio = 1.0
% xsize= 6.5
% ysize=xsize*aspectratio
% charsize = 1.0
% !p.charsize = charsize
% 
% PS_Start,filename=PlotFileName, /encapsulated, font=1, tt_font='Helvetica', $
%          charsize=charsize, default_thickness=2, /nomatch, xsize=xsize, ysize=ysize
%  
% ; pick color tables below, since may want different colors or ranges for different contour plots 
% 
% ; ***** COMPUTE THE 1-SIDED, 2D VARIANCE SPECTRUM Psi1S(kx1S,ky1S) ***** 

% ; Use the 1-sided, 2D variance spectrum using Psi(k,phi) of Elfouhaly et al Eq. (67) and 
% ; Psi(kx,ky) = k Psi(k,phi).  Psi1s is needed only for "downwind" to "crosswind" directions.

% print,' '
% print,format='(4x,a,f5.1)', 'The 2D variance spectrum is Elfouhaily et al. with U10 =',U10 
% print,' '
% ;----- Compute the 1-sided variance spectrum Psi1s(kx,ky) for all kx,ky values. 
% 
% ; The variance spectrum is computed only for "downwind" directions (kx1S,ky1S).  
% ; The upwind directions (negative kx) are obtained by symmetry: E(-kvec) = E(kvec).
% 
% ; To avoid having multiple large arrays, make multiple calls to ECKV2D_k_phi with 
% ; scalar k and phi values. This calculation is done only once if multiple random 
% ; surfaces are generated for a given variance spectrum.
% 
% ; Psi1s(kx,ky) = Psi(k,phi) for the corresponding (kx,ky) and (k,psi) points.
% 
% ; Note that Psi1s still has the magnitude of the one-sided variance spectrum as usually
% ; computed and plotted, even though it has been extended to negative frequencies.
% ; The factor of 1/2 on Psi1s to get a two-sided spectrum Psi2S is included in C3.

Psi1s=zeros(Nx,Ny);%holds all kxmath and kymath frequencies
for ikx=1:Nx/2 % loop over non-neg kx values kx1S 1-Nx/2, exclude 0 and go to kx1S(33)?
    for iky = Ny/2:Ny %non-negative ky values Ny/2 - Ny

          k = sqrt(kx1S(ikx).^2 + ky1S(iky).^2)

          phirad=atan2(ky1S(iky),kx1S(ikx));
          
          Psi1s(ikx+Nx/2,iky) = ECKV2D_k_phi(k,phirad,U10); % Psi1s(kx,ky) = Psi1s(k,phi)
          if iky >= Ny/2+1 && iky <= Ny-1
              Psi1s(ikx+Nx/2,Ny-iky+1)=Psi1s(ikx+Nx/2,iky);
          end
    end
end

%; fill in the negative kx values by symmetry
for iky=1:Ny
  for ikx = 1:Nx/2-1
      Psi1s(ikx,iky) = Psi1s(Nx-1-ikx,iky); 
  end
end
%; set (0,0) frequency to 0 to force MSL = 0:
Psi1s(Nx/2,Ny/2) = 0.0;


%; Psi1s now has all frequencies in math order
% ;print,' '
% ;print,' Psi1s(kxmath,kymath); (0,0) at center;  Psi1s(-kvec) = Psi1s(+kvec)'
% ;print,format='(a8,11x,16i12)',' ikx =',indgen(Nx)
% ;print,format='(a8,9x,16f12.3)',' kx =',kxmath
% ;for iky=Ny-1,0,-1 do begin
% ;    print,format='(a8,i5,f9.3,16e12.2)',' iky, ky =',iky,kymath(iky),Psi1s(*,iky)
% ;endfor
% 
% ; ***** CONTOUR THE ONE-SIDED, 2D VARIANCE SPECTRUM; UPPER LEFT PANEL *****  
% 
Psi1splot = Psi1s(Nx/2+1:Nx,:)'; %; temp array for plotting the right half of the symmetrical variance spectrum
figure(4)
subplot(2,2,1)
vpsi=[10^0, 10^-1, 10^-2, 10^-3 10^-4 10^-5 10^-6 10^-7 10^-8];
% [cpsi,hpsi]=contourf(((abs((Psi1splot)))),5);%vpsi)
% clabel(cpsi,hpsi,vpsi)
% ;Psi1splot = Psi1s  ; for contouring the 2sided spectrum (incl neg kxmath)
contourf(10*log10(Psi1splot))
colorbar
title('Spread Eckv elfouhaily')
subplot(2,2,2)

loglog(kx1S(1,1:32),Psi1splot(32,:))
title('Spectrum along the x-axis from surface')
xlabel('spatial frequency k [rad/m]')
ylabel('Variance spectrum S(k) [m^2/rad/m]')

subplot(2,2,3)
loglog(kx1S,elfouhailyD(U10,0.84,0,kx1S))
title('Spectrum from implemented elfouhaily')
xlabel('spatial frequency k [rad/m]')
ylabel('Variance spectrum S(k) [m^2/rad/m]')

subplot(2,2,4)
loglog(kx1S(1,1:32),Psi1splot(32,:))
hold on
loglog(kx1S,elfouhailyD(U10,0.84,0,kx1S))
title('Implementer + from surface')
xlabel('spatial frequency k [rad/m]')
ylabel('Variance spectrum S(k) [m^2/rad/m]')
legend('Andreas implementation','From surfacespectrum')

% surfc(1:64,1:64,real(Psi1s))

% ; ***** GENERATE A RANDOM HERMITIAN zhat(kx,ky) ***** 
% 
% ; Recall the note above on frequency order in FFTs.
%  
% ; When creating a random zhat for the IDL inverse FFT routine:
% ;   The order of the random amplitudes must be in the FFT frequency order;
% ;   The real and imag parts must give a Hermitian complex amplitude, 
% ;      so that the inverse FFT gives a real function
% ;   Hermitian means H(-k) = H*(k); H* = complex conjugate of H
% 
% ; The computation order is as follows:
%  
% ; (1) Draw random numbers and define a random, complex zhat(kx1S,ky1S) using the
% ;     the one-sided variance spectrum Psi1s(kx1S,ky1S) computed above.
%           
%           rho = normrnd(0,1,L,L);
%           sigma = normrnd(0,1,L,L);      
%
% ; (2) Use Hermitian symmetry to define zhat(kxmath,kymath) for the negative kxmath frequencies.
% ;     Note that variance spectra satisfy E(-kvec) = E(kvec) and E(kvec=0) = 0 ; kvec=(kxmath,kymath).
% ;     The math freq order is convenient for defining zhat and for contour plotting to check the symmetry of zhat.    
% ;     Check Hermitian requirements on Real{zhat} and Imag{zhat} by contouring these arrays:
% ;     Recall that Real{zhat} is an even functiion of kmath
% ;     Recall that Imag{zhat} is an odd functiion of kmath    
%     
%           zhat0 = (1./sqrt(2)).*(rho + j.*sigma).*sqrt(Psi1s).*exp(-j.*w.*t);
%               note that for t = 0 the exponential term reduces to 1.
%           zhat = (1./sqrt(2)).*(zhat0+ctranspose(zhat0));
%           Function ishermetian(zhat) = 1, means this is hermitian
%
% ; (3) Shift zhat(kxmath,kymath) to get zhat(kxFFT,kyFFT) for the FFT
% ;       is it supposed to be z = circshift(z,-Nx/2,-Ny/2) or
% ;       should we shift Psi1s to Psi2s? Lets try it out I guess...
% ; (4) compute z(x,y) = Real{invFFT{zhat(kxFFT,kyFFT)}}
% ;         Z = ifft(z);
% ;     Check that Imag{invFFT{zhat}} = 0
% ;         As of right now, imag(Z) != 0 but this might work it self out
% when we fixe the twosided spectrum with circshift.
% ;         
    

% ; shift from math to FFT frequency order

Psi1s = circshift( Psi1s,[Nx/2-1 Ny/2-1]); %Minus 0 or 2? Matlab note
                                                %should be N/2 + 1, thus
                                                %Psi1s(N/2) = 0 to get
                                                %"mathorder" so -N/2-1
                                                %should be good for FFT
                                                %order, right?

% ; (0,0) freq should now be at array location [0,0] OK as is right now
% ; print,'FFT freq order: Psi1s[0,0] (should = 0) = ',Psi1s[0,0]

% ; C3 is the product of
% ; the the missing 1/sqrt(2) from sqrt(Psi1s/2) in Tessendorf's Eq. 41
% ; the 1/sqrt(2) in the def of hhat_o in Tessendorf Eq 41
% ; the missing 1/sqrt(2) in Tessendorf's equation 43
C3 = 1/sqrt(8);

% ; Psi1s is the value of the one-sided, continuous variance spectrum at freq (kx,ky) in units of m^2/(rad/m)^2.  
% ; Must multiply by Deltakx and Deltakx to get variance (m^2) contained in the finite spatial 
% ; frequency intervals Deltakx and Deltaky.  See Light and Water Eq. (4.78).
% ; Psiroot is now the square root of the discrete TWO-sided spectrum, Psi_2s(u,v), including all factors 
% ; missing from the Tessnedorf equations.

Psiroot = C3*sqrt(Psi1s*Deltakx*Deltaky); 

% Psi1s = 0; %; now done with Psi1s array; free storage - Better to keep it
% for debugging purposes - Oskar


% ;***** CAN START LOOPING HERE TO GENERATE MULTIPLE SURFACES 
% ;      FOR THE SAME VARIANCE   SPECTRUM BUT DIFFERENT RANDOMIZATIONS.
% 
% 
% ; ----- Draw a pair of independent N(0,1) random numbers (ranr, rans) for each kx and ky
% 
% ; These will set the random amplitudes and phases of the waves, or real and imagnary parts of zhat
% 
% ; RANDOMU is the IDL random number generator (period ~ 10^8); /normal gives normal(0,1) random variables
ranr=zeros(Nx,Ny); %; for the random real parts of zhat
rans=zeros(Nx,Ny); %for the random imaginary parts 

for ikx=1:Nx
    ranr(ikx,:)=randn(1,Ny);
end
for ikx=1:Nx
    rans(ikx,:)=randn(1,Ny);
end

%; define the complex Hermitian zhat array for the full range of pos and negative frequencies
zhat=ones(Nx,Ny);

% ;----- Compute zhat(kx,ky)
% 
% ; NOTE: this code computes zhat at time = 0, which gives cos(omega t) = 1 and sin(omega t) = 0
% ; in the general time-dependent equations used in cgAnimate2D.pro.  The zhat equation below is
% ; therefore not the most general form for zhat(kx,ky,t).
% 
% ; The general method is based on Tessendorf (2004; Simulating Water Surfaces), as 
% ; corrected for consistency with the chosen one-sided wave variance spectrum
% 
% ; note that for the FFT frequency ordering, k(-j) = k(N-j)
% 
% ; Must now loop over all frequencies, and treat the Nyquist frequencies as special cases
% 
% ; the non-zero and non-Nyquist frequencies:
% Note from Andreas: MATLAB work from 1->Nx unlike IDL 0->Nx-1

%Unsure on the indexing - Andreas 
for ikx=2:Nx/2
    for iky=2:Ny
        zhat(ikx,iky) =(ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky))+ j*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky));
      zhat(Nx+2-ikx,Ny+2-iky) = conj( zhat(ikx,iky) );
    end
end


% ; all ky for kx = 0 freq at index 0 and the Nyquist frequency at kx index Nx/2:
 for iky=2:Ny/2
    ikx = Nx/2+1;
    zhat(ikx,iky) = (ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky))+i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky));
    zhat(ikx,Ny+2-iky) = conj( zhat(ikx,iky) );
    ikx = 1;
    zhat(ikx,iky) = (ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(ikx,Ny+2-iky) * Psiroot(ikx,Ny+2-iky))+i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(ikx,Ny+2-iky) * Psiroot(ikx,Ny+2-iky));
    zhat(ikx,Ny+2-iky) = conj( zhat(ikx,iky));
 end
%  Temp sol for NaN in Ps1s
%  zhat(2,3)=7*10^-6;
%  zhat(5,3)=7*10^-6;
%  zhat(Nx,Ny-1)=7*10^-6;
%  zhat(Nx+2-5,Ny+2-3)=7*10^-6;
%  ; all kx for ky = 0 and Nyquist frequency at ky index Ny/2:
 for ikx=2:Nx/2+1
   iky = Ny/2+1;
    zhat(ikx,iky) =(ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky))+i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky));
    zhat(Nx+2-ikx,iky) = conj( zhat(ikx,iky) );
   iky = 1;
    zhat(ikx,iky) = ( ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,iky) * Psiroot(Nx+2-ikx,iky))+i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,iky) * Psiroot(Nx+2-ikx,iky));
    zhat(Nx+2-ikx,iky) = conj( zhat(ikx,iky) );
 end
 
%  ; Nyquist ky freq at ky index Ny/2
 ikx = 1;
 iky = Ny/2+1;
  zhat(ikx,iky) = ( ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(ikx,Ny+2-iky) * Psiroot(ikx,Ny+2-iky))+i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(ikx,Ny+2-iky) * Psiroot(ikx,Ny+2-iky));
% ; Nyquist kx freq at kx index Nx/2
 ikx = Nx/2+1;
 iky = 1;
 zhat(ikx,iky) = ( ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,iky) * Psiroot(Nx+2-ikx,iky)) +i*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,iky) * Psiroot(Nx+2-ikx,iky));

% ; the "double" Nyquist frequency at Nx/2,Ny/2
 ikx = Nx/2+1;
 iky = Ny/2+1;
 zhat(ikx,iky) = ( ranr(ikx,iky) * Psiroot(ikx,iky) + ranr(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky))+ j*(rans(ikx,iky) * Psiroot(ikx,iky) - rans(Nx+2-ikx,Ny+2-iky) * Psiroot(Nx+2-ikx,Ny+2-iky));
 
%  ; set the (0,0) value to 0 (MSL = 0)
 zhat(1,1) =(0+ i*0);
 
%  ; zhat as just defined has the frequencies in FFT order, as needed by the FFT routine.
% 
% ; for plotting the real and imaginary parts in math order, shift back to math order:
% ; The math order can be obtained from the FFT order by a circular shift by N/2-1 to the right:
% ; kmath = SHIFT(kFFT,Nx/2-1)

Zreal = circshift( real(zhat), [Nx/2, Ny/2]);
Zimag = circshift( imag(zhat), [Nx/2, Ny/2]);
idebug = 0;

if idebug == 1 
%{
print,' '
print,' Real{zhat}(kxmath,kymath); (0,0) at center;  Real{zhat(-k)} = Real{zhat(k)}'
print,format='(a8,11x,16i12)',' ikx =',indgen(Nx)
print,format='(a8,9x,16f12.3)',' kx =',kxmath
for iky=Ny,1
    print,format='(a8,i5,f9.3,16e12.2)',' iky, ky =',iky,kymath(iky),zreal(*,iky)
end

print,' '
print,' Imag{zhat}(kxmath,kymath); (0,0) at center;  Imag{zhat(-k)} = -Imag{zhat(k)}'
print,format='(a8,11x,16i12)',' ikx =',indgen(Nx)
print,format='(a8,9x,16f12.3)',' kx =',kxmath
for iky=Ny-1,0,-1 do begin
    print,format='(a8,i5,f9.3,16e12.2)',' iky, ky =',iky,kymath(iky),zimag(*,iky)
endfor

%}
end

vzplot=[0.02 0.016 0.012 0.008 0.004 0 -0.004 -0.008 -0.012 -0.016 -0.020];
figure(5)
contourf(linspace(-2,2,64),linspace(-2,2,64),Zreal',vzplot)
colorbar

title('zhat real')

figure(6)
contour(Zimag')
title('zhat real no shift')
contourf(linspace(-2,2,64),linspace(-2,2,64),real(zhat)',vzplot)


% ; ***** TAKE THE INVERSE FFT TO GET THE SEA SURFACE *****
disp('Frequency-shifted zhat(0,0) [should = (0,0)] ='+num2str(zhat(1,1)));

% ; take the inverse FFT of zhat in FFT frequency order

zcomplx = ifft2(zhat);

% ; ***** EXTRACT THE SEA SURFACE *****
% 
% ; The surface wave ampltudes are the real part of the inverse FFT

zsurf = real(zcomplx);
figure(7)
vzsurf=[0.6 0.45 0.3 0.15 0 -0.15 -0.3 -0.45 -0.6]
ZSURF = 64*64*zsurf;
surfc(linspace(0,100,64),linspace(0,100,64),ZSURF)
colorbar
zimag = imag(zcomplx);

% ; ----- Checks on the generated surface
disp('Checks on the generated z(x,y):')

% ; compute the mean sea surface, which should be zero, the mean squared surface elevation,
% ; and sum of amplitudes squared from full (two-sided) zhat array
sumz = 0;
sumzsq = 0;
sumzhat2 = 0;
for ix = 1:Nx
  for iy = 1:Ny
    sumz = sumz + zsurf(ix,iy);
    sumzsq = sumzsq + zsurf(ix,iy).^2;
    sumzhat2 = sumzhat2 + abs(zhat(ix,iy))^2; %; (ix,iy) ranges same as (u,v) ranges
  end
end

% ; compute total energy from the one-sided (1S) energy spectrum
% ; add terms that have only 0 or pos freqs once, and double terms than have both pos and neg freqs
sumzhatsq = abs(zhat(1,1)).^2;
for iky=1:Ny  
    sumzhatsq = sumzhatsq + abs(zhat(1,iky)).^2 + abs(zhat(Nx/2+1,iky)).^2; %; zero and Nyquist in x for all y; occur once
end
 for ikx=2:Nx/2
    for iky=1:Ny
    sumzhatsq = sumzhatsq + 2.0*abs(zhat(ikx,iky))^2; %; all other freqs occur twice as pos and neg freqs
    end
 end


 zavg = sumz/(Nx*Ny);

%  ; FIX THIS--BAD VALUES
% disp('   avg z(x,y) (should = 0) = '+num2str(zavg));
% disp('   max |Imag{zcomplx}| (should = 0) ='+num2str(max(abs(imag(zcomplx)))));
% 
% disp('Parsevals identity:')
% disp('                            sum z^2 ='+num2str(sumzsq));
% disp('   Nx * Nx * sum zhat^2 (one-sided) ='+num2str(Nx*(Ny)*sumzhatsq));
% disp('   Nx * Nx * sum zhat^2 (two-sided) ='+num2str((Nx)*(Ny)*sumzhat2));


% ; the significant wave height from avg of z^2
avgzsq = sumzsq/((Nx)*(Ny));
H13 = 4.0*sqrt(avgzsq); %; approximate significant wave height
% disp('   significant wave height H_1/3 = 4 SQRT(avg{z^2}) ='+num2str(H13));

% ; compute the mean square slopes and mean slope angle in x and y directions
dzdx2 = 0;
thetax = 0;
for ix = 2:Nx
  for iy = 1:Ny
    dzdx2 = dzdx2 + (zsurf(ix,iy) - zsurf(ix-1,iy))^2;
    thetax = thetax + abs(atan2(zsurf(ix,iy) - zsurf(ix-1,iy),Deltax)); %removed radeg matlab inputs in rad for atan default
  end
end
dzdx2 = dzdx2/(Deltax*Deltax*(Nx)*(Ny+1));
thetax = thetax/((Nx)*(Ny+1));
%{
print,'   sample mean square slope, alongwind = ', dzdx2
print,'   sample mean square slope, crosswind = ', dzdy2
print,'   total sample mean square slope      = ',dzdx2 + dzdy2
print,'   Cox-Munk mss_x = ',0.0316*U10
print,'   Cox-Munk mss_y = ',0.0192*U10
print,'   Cox-Munk mss   = ',0.001*(3.0 + 5.12*U10)
print,'   sample avg slope angle, alongwind = ', thetax
print,'   sample avg slope angle, crosswind = ', thetay
%}

% zcomplx = 0 %; now done with zcomplx array; free storage

%  Save surface to .mat file
SurfaceSave=zeros(length(Zreal(:,1))*length(Zreal(1,:)),3);
surlen=length(Zreal(1,:));

%for a=1:surlen
%   for b=1:surlen
%        
%   end
%end
    
    
for a=1:length(Zreal(:,1))
        for c=1:length(Zreal(:,1))
            SurfaceSave((a-1)*length(Zreal(:,1))+c,1)=c*Deltax-50;
            SurfaceSave((a-1)*length(Zreal(:,1))+c,2)=a*Deltay-50;
            SurfaceSave((a-1)*length(Zreal(:,1))+c,3)=Zreal(c,b);
        end
end
        
disp('Save surface')

save('surface_elfouhaily.mat','SurfaceSave');
save_check=load('surface_elfouhaily.mat','SurfaceSave')

% save('surface_elfouhaily.mat','SurfaceSave');
%  save('MyMatrix.txt', 'SurfaceSave', '-ascii', '-float', ' ')

% save_check=load('surface_elfouhaily.mat','SurfaceSave')

 dlmwrite('myFile.txt',SurfaceSave,'delimiter',' ');
