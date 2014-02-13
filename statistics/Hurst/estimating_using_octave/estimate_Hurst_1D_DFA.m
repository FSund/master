function H = estimate_Hurst_1D_DFA(X)

    [N, M] = size(X);
    if (M > 1 && N > 1)
        error('Input needs to be a 1d vector. Aborting!')
    end
    if (M > 1)
        X = X'; % transpose to get column vector in case input is row vector
    end
    N = max(M,N);
    smin = 6;               % approx, can be tuned
    smax = floor(N/4.0);    % approx, can be tuned
%     svec = smin:smax;
    svec = round(logspace(log10(smin), log10(smax)));
    svec = unique(svec);    % remove duplicates
    
    F2 = zeros(1, length(svec));
    for is = 1:length(svec) % loop over sizes
        s = svec(is);
        Ns = floor(N/s); % number of segments
        F2sum = 0.0;
        for v = 1:Ns % loop over segments
            ivec = (v-1)*s + (1:s);

            % Cumulative sum
            uvw = cumsum(X(ivec));

            % Fit to different polynomials
            x = ivec';  % transpose to get column vector
            y = uvw;    % already column vector

            % Select polynomial
            % % eq (2), ai + bj +c
            A = [x ones(size(x))];
            % % eq (3), ai^2 + bj^2 + c
%             A = [x.*x ones(size(x))];

            % Find coefficients using solve(A,B) == A\B
            coeff = A\y;
            U = A*coeff;

            % Find residual matrix eps_{u,v}
            eps = uvw - U;

            % Detrended fluctuation function F^2(u,w,s), eq (8)
            F2sum = F2sum + sum(eps.*eps)/s;
        end
        % Overall detrended fluctuation, eq. (9)
        F2(is) = F2sum/Ns;
    end
    
    x = log10(svec);
    y = log10(sqrt(F2));
    p = polyfit(x,y,1);
    H = p(1);
%     fprintf('H = %1.4f\n', H);
% 
%     figure;
%     plot(x, y);
%     hold all;
%     x = log10([smin/2 smax*2]);
%     y = polyval(p, x);
%     plot(x, y);
%     legend('log10(sqrt(FF))', sprintf('%1.4fx + %1.4f', p(1), p(2)), 'Location', 'Best');
end
