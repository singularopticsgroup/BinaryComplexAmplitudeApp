function [LG]  = Laguerre_Gaussian(Wavelength,p,l,w0,pixel_size)

% Wavelength 
% p - radial index
% l - azimuthal index
% w0 - beam waist
% pixel_size - size of DMD pixel

W = 1920; % Width
H = 1080; % Height
x = -W/2:1:(W/2-1); 
y = -H/2:1:(H/2-1);
px = pixel_size*10^(-6);
x1 = x.*px;
y1 = y.*px;
    
% Creates matrices with coordination system
[X1, Y1] = meshgrid(x1,y1);

                [phi, r] = cart2pol(X1, Y1);
                k = 2*pi/Wavelength;                  % Wavenumber of light
                zR = k*w0^2/2;
                z=0;

                 w = w0 * sqrt(1 + z.^2/zR^2);
                 R = sqrt(2)*r./w;
                % 
                % Lpl from OT toolbox 
                Lpl = nchoosek(p+l,p) * ones(size(R));   % x = R(r, z).^2

                for t = 1:p
                    Lpl = Lpl + (-1)^t/factorial(t) * nchoosek(p+l,p-t) * R.^(2*t);
                end

                LGN = @(m,p)1/sqrt(2)*(r/w0).^abs(m).*exp(-r.^2/w0.^2).*Lpl.*exp(1i*m*phi).*exp(-1i*(2*p + m + 1)*atan(z/zR));
                LG = LGN(l,p);
end