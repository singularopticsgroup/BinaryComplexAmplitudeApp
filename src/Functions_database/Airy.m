function A = Airy(Wavelength,b,s,pixel_size)            

% b - phase factor
% s - scale
% pixel_size - size of DMD pixel

W = 1920;     % Width
H = 1080;     % Height
x = -W/2:1:(W/2-1); 
y = -H/2:1:(H/2-1);
px = pixel_size*10^(-6);
x1 = x.*px;
y1 = y.*px;
    
% Creates matrices with coordination system
[X1, Y1] = meshgrid(x1,y1);

            l = Wavelength*1e-9;
            n = 1.2;
            alpha = 0.2;
            alpha=alpha*pi/180;
                                
            k=2*pi/l;
            
            kr=k*alpha*(n-1);
            A = 4*airy(X1.*kr/s).*airy(Y1.*kr/s).*exp(b*(X1/s+Y1/s).*kr);

end