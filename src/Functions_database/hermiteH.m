function H_p = hermiteH(n, x)

                % generate the nth hermite polynomial of x

                
                m1 = 0:floor(n/2);
                % s = numel(x);
                [g1, g2] = size(x);
                %X=repmat(x,[ones(1,s),g1, g2]);
                Xh = repmat(x, [1, numel(m1)]);
                m1 = repmat(m1, [g1, g2]);
            
                H_p = factorial(n) * sum((-1).^m1 ./ (factorial(m1) .* factorial(n - 2*m1)) .* (2*Xh).^(n - 2*m1), 2);
end