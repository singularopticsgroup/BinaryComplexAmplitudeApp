function zernike_map = do_zernike(w,h,n,m,coefficient)
            
%Functon will calculate proper zernike polynomial Z_n_m and multiply it by specified coeffcient
% w,h - resolution (x,y)

            if coefficient~=0 %In case of 0 value, returns an empty map  

                x = linspace(-1,1,h);
                [X_z,Y_z] = meshgrid(x,x);
                [theta,r] = cart2pol(X_z,Y_z);

                 idx = r<=1;
                 z = nan(size(X_z));
                
                % calculate polynomial
                z(idx) = zernike_function(m,n,r(idx),theta(idx));
                z = padarray(z,[0 (w-h)/2], 0, 'both');
                z(isnan(z))=0;
                zernike_map=mod(z*coefficient,2*pi);
            else
                zernike_map=zeros(h,w);
            end

end