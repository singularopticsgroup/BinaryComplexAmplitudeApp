% EVEN Ince Polynomial

function [IP,eta,coef,dIP]=CInceIGB(p,m,q,z) 
            %{
            This function calculates the EVEN Ince polynomials which are solutions of the Ince equation.
            The Ince equation is a periodic linear second-order differential equation that has two families of independent solutions, 
            namely, the even C^{m}_p(z,q) and odd S^{m}_p(z,q) Ince polynomials of order p and degree m. 
            Physical considerations are such that Ince polynomials are periodic with period 2p. The values of \etha that satisfy this condition
            are the eigenvalues of the Ince equation which is given by
            
            d^2F/dz^2 + q*sin(2z)dF/dz+(\eta-p*q*cos(2z))*F=0
            
            where
            0<=z<2pi
            p=0,1,2,3,...  
            q=complex parameter
            \etha = eigenvalue of the Ince equation.
            
            INPUTS:
                p=0,1,2,3... Order of Ince polynomial
                0<=m<=p      m is the degree of the Ince polynomial
                (p,m) must have the same parity, i.e., (-1)^(p-m)=1
                q = complex parameter
                0<=z<2pi     independent variable (Vector or Matrix)
            
            OUTPUTS:
                IP=C^{m}_p(z,q) EVEN Ince Polynomial
                eta  = eigenvalue of the Ince Polynomial 
                coef = coefficients of the Ince Polynomial
                dIP  = derivative of the Ince Polynomial
            %}
            
            % Check Input
            if (m<0)||(m>p); error('ERROR: Wrong range for "m", 0<=m<=p'); end
            if (-1)^(m-p)~=1;   error('ERROR: (p,m) must have the same parity, i.e., (-1)^(m-p)=1');  end
            [largo,ancho]=size(z); % change input to vector format
            z=transpose(z(:));
            normalization=1;
            
            % Calculate the Coefficients 
            if mod(p,2)==0
                %%%% p Even %%%%
                j=p/2;  N=j+1;  n=m/2+1; 
                
                % Matrix
                M=diag(q*(j+(1:N-1)),1) + diag([2*q*j,q*(j-(1:N-2))],-1) + diag([0,4*((0:N-2)+1).^2]);
                if p==0; M=0; end
                    
                % Eigenvalues and Eigenvectors 
                [A,ets]=eig(M);
                ets=diag(ets); 
                [ets,index]=sort(ets);
                A=A(:,index);
                
                % Normalization
                if normalization==0  
                   N2=2*A(1,n).^2+sum(A(2:N,n).^2);
                   NS=sign(sum(A(:,n)));
                   A=A/sqrt(N2)*NS;
                else 
                   mv=(2:2:p).';
                   N2=sqrt(A(1,n)^2*2*gamma(p/2+1)^2+sum((sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(2:p/2+1,n)).^2 ));
                   NS=sign(sum(A(:,n)));
                   A=A/N2*NS;
                end
                
                % Ince Polynomial
                r=0:N-1;
                [R,X2]=meshgrid(r,z);
                IP=cos(2*X2.*R)*A(:,n);
                dIP=-2*R.*sin(2*X2.*R)*A(:,n);
                eta=ets(n);
                    
            else
                %%%% p ODD %%%
                j=(p-1)/2;  N=j+1;  n=(m+1)/2; 
                
                % Matrix
                M=diag(q/2*(p+(2*(0:N-2)+3)),1)+diag(q/2*(p-(2*(1:N-1)-1)),-1) + diag([q/2+p*q/2+1,(2*(1:N-1)+1).^2]);
            
                % Eigenvalues and Eigenvectors 
                [A,ets]=eig(M);
                ets=diag(ets);
                [ets,index]=sort(ets);
                A=A(:,index);
            
                % Normalization
                if normalization==0  
                   N2=sum(A(:,n).^2);
                   NS=sign(sum(A(:,n)));
                   A=A/sqrt(N2)*NS;
                else
                    mv=(1:2:p).';
                    N2=sqrt(sum( ( sqrt(gamma((p+mv)/2+1).*gamma((p-mv)/2+1) ).*A(:,n)).^2 ));
                    NS=sign(sum(A(:,n)));
                    A=A/N2*NS;
                end
                
                % Ince Polynomial
                r=2*(0:N-1)+1;  
                [R,X2]=meshgrid(r,z);
                IP=cos(X2.*R)*A(:,n);
                dIP=-R.*sin(X2.*R)*A(:,n);
                eta=ets(n);
            end
            
            coef=A(:,n);
            IP=reshape(IP,[largo,ancho]); % reshape output to original format
            dIP=reshape(dIP,[largo,ancho]);
        
end