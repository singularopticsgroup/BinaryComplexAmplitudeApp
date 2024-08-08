function [IGB,X,Y]=Ince_Gaussian(L,N,parity,p,m,e,w0,k,z)
            
% The mesh_elliptic, CInceIGB, SInceIGB must be in the same folder as Ince_Gaussian function
            %{
            Calculate an Ince-Gaussian Beam at a given z plane
            
            INPUTS:
                L = transverse physical size of the X-Y space [-L,L,-L,L]
                N = number of sampling points (must be ODD)
                parity = parity of the beam, 0 = EVEN  C;  1 = ODD  S
                p, m = order and degree of the Ince Gaussian beam
                       p=0,1,2,3,... for parity=0 and p=1,2,3,.. for parity=1
                       0<=m<=p   
                       (p,m) must have the same parity, i.e., (-1)^(p-m)=1
                e = ellipticity parameter
                w0 = beam width(waist) at z=0
                k = 2*pi/lambda, wavenumber, lambda=wavelength 
                z = propagation distance
            
            OUTPUTS:
                IGB = Ince Gaussian Beam (Intensity is analytically normalized to 1, i.e., Integral(|IGB|^2)=1)
                X and Y = space matrices
                
            %}
            
            % CHECK INPUT
            if mod(N,2)==0;      error('ERROR: N must be ODD'); end
            if parity==0
                if (m<0)||(m>p); error('ERROR: Wrong range for "m", 0<=m<=p'); end
            else
                if (m<1)||(m>p); error('ERROR: Wrong range for "m", 1<=m<=p'); end  
            end
            if (-1)^(m-p)~=1;    error('ERROR: (p,m) must have the same parity, i.e., (-1)^(m-p)=1');  end
            
            % PARAMETERS
            f0 = sqrt(e/2)*w0; % Focal distance of elliptic coordinates at z=0
            
            % INCE GAUSSIAN BEAM
            if z==0 % At the z=0
               [xhi,etha,X,Y]=mesh_elliptic(f0,L,N);% Calculate elliptic coordinates (xhi,etha), (X,Y)= Cartesian Coordinates
               R = sqrt(X.^2 + Y.^2);                   % Radial coordinate
               if parity==0
                  IGB=CInceIGB(p,m,e,etha).*CInceIGB(p,m,e,1i*xhi).*exp(-(R/w0).^2);  % Even IGB
               else
                  IGB=SInceIGB(p,m,e,etha).*SInceIGB(p,m,e,1i*xhi).*exp(-(R/w0).^2);  % Odd IGB
               end
            else % At z~=0
               zr=1/2*k*w0^2;                           % Rayleigh range 
               wz=w0*sqrt(1+(z/zr).^2);                 % Beam width(waist) at z=0   
               Rz=z*(1+(zr./z).^2);                     % Radius of curvature of the phase front
               f = f0*wz/w0;                            % Focal distance of elliptic coordinates at z   
               [xhi,etha,X,Y]=mesh_elliptic(f,L,N); % Calculate elliptic coordinates (xhi,etha), (X,Y)= Cartesian Coordinates
               R = sqrt(X.^2 + Y.^2);                   % Radial coordinate
               if parity==0
                  IGB=(w0/wz)*(CInceIGB(p,m,e,etha).*CInceIGB(p,m,e,1i*xhi)).*exp(-(R/wz).^2).* ...
                     exp(1i*(k*z + k*R.^2/(2*Rz)-(p+1)*atan(z/zr)));    % Even IGB
               else
                  IGB=(w0/wz)*(SInceIGB(p,m,e,etha).*SInceIGB(p,m,e,1i*xhi)).*exp(-(R/wz).^2).* ...
                     exp(1i*(k*z + k*R.^2/(2*Rz)-(p+1)*atan(z/zr)));    % Odd IGB
               end   
            end
            
            % COMPUTE THE NORMALIZATION CONSTANTS
            if parity == 0  
               if mod(p,2)==0   
                  [C0,~,coef]=CInceIGB(p,m,e,0);
                  [Cp,~,~]=CInceIGB(p,m,e,pi/2);
                  Norm = (-1)^(m/2)*sqrt(2)*gamma(p/2+1)*coef(1) *sqrt(2/pi)/w0/C0/Cp;  % Calculate normalization constant
               else
                  [C0,~,coef]=CInceIGB(p,m,e,0); 
                  [~,~,~,DCp]=CInceIGB(p,m,e,pi/2);
                  Norm = (-1)^((m+1)/2) * gamma((p+1)/2+1) * sqrt(4*e/pi) * coef(1) / w0 / C0 / DCp; % Calculate normalization constant
               end
            else             
               if mod(p,2)==0   
                  [~,~,coef,dS0]=SInceIGB(p,m,e,0);
                  [~,~,~,dSp]=SInceIGB(p,m,e,pi/2); 
                  Norm = (-1)^(m/2)*sqrt(2)*e*gamma((p+2)/2+1)*coef(1) *sqrt(2/pi)/w0/ dS0 / dSp;  % Calculate normalization constant      
               else
                  [Sp,~,coef,~]=SInceIGB(p,m,e,pi/2); 
                  [~,~,~,dS0]=SInceIGB(p,m,e,0);       
                  Norm = (-1)^((m-1)/2) * gamma((p+1)/2+1) * sqrt(4*e/pi) * coef(1) / w0 / Sp / dS0;  % Calculate normalization constant      
               end
            end
            
            IGB = IGB*Norm; % Normalized the IGB to Integral(|IGB|^2)=1
        
end