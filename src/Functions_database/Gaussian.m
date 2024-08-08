function G = Gaussian(pixel_size,w0)

% w0 - beam waist
% pixel_size - size of DMD pixel

W = 1920;   % Width
H = 1080;   % Height
x = -W/2:1:(W/2-1); 
y = -H/2:1:(H/2-1);

% Creates matrices with coordination system
[X, Y] = meshgrid(x,y);

            Xb = X*pixel_size*1e-6;
            Yb = Y*pixel_size*1e-6;
            rho = sqrt(Xb.^2+Yb.^2);

            G = exp(-rho.^2/w0^2);
end