function [B] = Bessel(l,n,alpha,m1,pixel_size)
                
% Bessel beam

%   l - Wavelength [m] 
%   n - Refractive index
%   alpha - Axicon angle [deg]
%   m1 - topological charge
% pixel_size - size of DMD pixel

W = 1920;   % Width
H = 1080;   % Height
x = -W/2:1:(W/2-1); 
y = -H/2:1:(H/2-1);
px = pixel_size*10^(-6);
x1 = x.*px;
y1 = y.*px;
    
% Creates matrices with coordination system
[X1, Y1] = meshgrid(x1,y1);

                l = l*1e-9;
                alpha = alpha*pi/180;

                Z=X1+1i*Y1;
    
                r=sqrt((X1.^2 + Y1.^2));         % Radius
                phi=angle(Z);         

                k=2*pi/l;

                kr=k*alpha*(n-1);
                kz=1;
                z=0;

                Bess=@(m1)exp(1i.*kz.*z).*besselj(m1,kr.*r).*exp(1i.*m1.*phi); 

                B = Bess(m1);
                
end