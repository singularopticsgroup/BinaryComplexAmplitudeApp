clear all; clc;
Wavelength = 633;           % Wavelength
p = 1;                      % p - radial index
l = 4;                      % l - azimuthal index
pixel_size = 8;             % pixel_size - size of DMD pixel
n = 1.2;                    % n - Refractive index
alpha = 0.1;                % alpha - Axicon angle [deg]
m1 = 2;                     % Topological charge for Bessel beam
L = 7.5e-3;                 % L = transverse physical size of the X-Y space [-L,L,-L,L]
N = 1081;                   % N = number of sampling points (must be ODD)
parity = 0;                 % Parity of Ince-Gaussian beam
p = 4;                      % angular eliptic mode
m = 2;                      % radial eliptic mode
e = 2;                      % ellipticity parameter
w0 = 1*10^(-3);             % beam waist
k = 2*pi/(Wavelength*1e-9); % wave number
z = 0;                      % z = propagation distance
Gi = 2;                     % Order of Super-Gaussian beam
Hi = 1080;                  % Height (Y resolution) for Hermite-Gaussian beam

% LG = Laguerre_Gaussian(Wavelength,p,l,w0,pixel_size);
% LC = Laguerre_Cosinus(Wavelength,p,l,w0,pixel_size);
% LS = Laguerre_Sinus(Wavelength,p,l,w0,pixel_size);
% B = Bessel(Wavelength,n,alpha,m1,pixel_size);
% [IGB,X,Y]=Ince_Gaussian(L,N,parity,p,m,e,w0,k,z);
% P = zeros(1080,1920);
% C=abs(IGB/max(max(IGB)));
% P(1:1080,(1920-1080)/2:(1920-1080)/2+1080-1) = C(1:1080,1:1080);

%G = Gaussian(pixel_size,w0);
%SG = Super_Gaussian(pixel_size,w0,Gi);
H_G = Hermite_Gaussian(pixel_size,0,1,Hi,w0);
HG = zeros(1080,1920);
HG(1:1080,(1920-1080)/2:(1920-1080)/2+1080-1) = H_G(1:1080,1:1080);

%A = Airy(633,0.1,3,8);

% ZERNIKE
%Functon will calculate proper zernike polynomial Z_n_m and multiply it by specified coeffcient
% w,h - resolution (x,y)
% n,m - indices eg. [0, 2] for defocus

%zernike_map = do_zernike(w,h,n,m,coefficient);
H = binarization(HG, 100, 100);
imshow(H)






