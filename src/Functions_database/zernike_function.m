% ZERNIKE Polynomials

function z = zernike_function(n,m,r,theta)
             
            n = n(:);
            m = m(:);
            r = r(:);
            theta = theta(:);
            length_r = length(r);
            m_abs = abs(m);
            rpowers = [];
            
            for j = 1:length(n)
                rpowers = [rpowers m_abs(j):2:n(j)];
            end
            rpowers = unique(rpowers);

            if rpowers(1)==0
                rpowern = arrayfun(@(p)r.^p,rpowers(2:end),'UniformOutput',false);
                rpowern = cat(2,rpowern{:});
                rpowern = [ones(length_r,1) rpowern];
            else
                rpowern = arrayfun(@(p)r.^p,rpowers,'UniformOutput',false);
                rpowern = cat(2,rpowern{:});
            end
            
            % Polynomials:

            z = zeros(length_r,length(n));
            for j = 1:length(n)
                s = 0:(n(j)-m_abs(j))/2;
                pows = n(j):-2:m_abs(j);
                for k = length(s):-1:1
                    p = (1-2*mod(s(k),2))* ...
                               prod(2:(n(j)-s(k)))/              ...
                               prod(2:s(k))/                     ...
                               prod(2:((n(j)-m_abs(j))/2-s(k)))/ ...
                               prod(2:((n(j)+m_abs(j))/2-s(k)));
                    idx = (pows(k)==rpowers);
                    z(:,j) = z(:,j) + p*rpowern(:,idx);
                end     
            end
            
            % Zernike functions:

            idx_pos = m>0;
            idx_neg = m<0;
            if any(idx_pos)
                z(:,idx_pos) = z(:,idx_pos).*cos(theta*m_abs(idx_pos)');
            end
            if any(idx_neg)
                z(:,idx_neg) = z(:,idx_neg).*sin(theta*m_abs(idx_neg)');

            end
end