function H_G = Hermite_Gaussian(pixel_size,n,m,Hi,w0)

% H - Height (Y resolution)
% n,m  horizontal and vertical Hermite-Gussian mode index
% pixel_size - size of DMD pixel
% w0 - beam waist

                %x = -Width.Value/2:1:(Width.Value/2-1); 
                y = -Hi/2:1:(Hi/2-1);
                %x1 = x.*pixel_size*10^(-6);
                y1 = y.*pixel_size*10^(-6);
                
                % Creates matrices with coordination system
                [Xh, Yh] = meshgrid(y1,y1);

                r = sqrt(Xh.^2+Yh.^2);
    
                E0 = exp(-r.^2 / w0^2); % create the gaussian
                
                H = hermiteH(n, sqrt(2) * Yh / w0)' .* hermiteH(m, sqrt(2) * Yh / w0) .* E0;

                H_G = abs(H);
                H_G = H_G ./ max(max(H_G));
            
end
