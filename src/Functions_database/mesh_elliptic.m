% Elliptic Coordinates
        function [xi,eta,X,Y]=mesh_elliptic(f,L,N)
            %{
            Creates an elliptic coordinate mesh of size 2Lx2L with N points
            
            INPUT:
                f = focal distance of elliptical coordinates 
                L = spatial size of the mesh, i.e., [-L,L]x[-L,L]
                N = ODD number of discrete points in the mesh. Must be ODD.
            
            OUPUT:
               xi  = elliptic coordinate mesh 
               eta = eta elliptic coordinate mesh
               X   = x-Cartesian coordinate mesh
               Y   = y-Cartesian coordinate mesh
            %}
            
            % Cartesian Coordinates
            [X,Y]=meshgrid(linspace(-L,L,N));
            
            % Elliptic Coordinates
            xi=zeros(N); eta=zeros(N);  % Initialization
            
            % Calculate First Quadrant             
            en = acosh((X(1:(N+1)/2,(N+1)/2:N)+1i*Y(1:(N+1)/2,(N+1)/2:N))/f);
            ee = real(en);
            nn = imag(en);
            nn = nn + (nn<0)*2*pi; 
            xi(1:(N+1)/2,(N+1)/2:N)=ee;
            eta(1:(N+1)/2,(N+1)/2:N)=nn;
                  
            % Calculate other quadrants by symmetry       
            xi(1:(N+1)/2 , 1:(N-1)/2)=fliplr(xi(1:(N+1)/2,(N+3)/2:N));
            xi((N+3)/2:N,1:N)=flipud(xi(1:(N-1)/2,1:N));
            eta(1:(N+1)/2,1:(N-1)/2)=pi-fliplr(eta(1:(N+1)/2,(N+3)/2:N));
            eta((N+3)/2:N,1:N)=pi+rot90((eta(1:(N-1)/2,1:N)),2);
        
        end