function H = binarization(Ew,nx,ny)
    %Ew - optical field, complex amplitude
    %nx - no. of grooves in x direction
    %ny - no. of grooves in y direction
    W = 1920;    % Width
    H = 1080;    % Height
    x = -W/2:1:(W/2-1); 
    y = -H/2:1:(H/2-1);
    
    % Creates matrices with coordination system
    [X, Y] = meshgrid(x,y);
    
    grating_phase = Grating(nx,ny,W,H,Y,X);

    % Lee
    b = asin(abs(Ew)/max(abs(Ew),[],'all'));

    H = 0.5+0.5.*(b+abs(Ew)).*sign(cos(angle(Ew)+grating_phase)-cos(b));
    H(H>0.5)=1;
    H(H<=0.5)=0;

end