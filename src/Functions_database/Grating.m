function [GP] = Grating(nx,ny,W,H,G_Y,G_X)

                % nx     Number of grooves in x direction
                % ny     Number of grooves in y direction
                % W      Width
                % H      Height
                % G_Y    Grating Y coordinate 
                % G_X    Grating X coordinate
 
                gx=nx/W; gy=ny/H;                            % Spatial freq
                
                GP = mod(2*pi*(G_Y*gy+G_X*gx),2*pi);         % Grating

                %GP=abs((angle(exp(1i*sin(2*pi*G_Y.*gy+2*pi*G_X.*gx)))+1)*pi-0.00001);

end