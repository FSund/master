function H = estimate_Hurst_1D_DFA(X)

    [N, M] = size(X);
    if (M > 1 && N > 1)
        error('Input needs to be a 1d vector. Aborting!')
    end
    N = max(M,N);
    smin = 6;               % approx, can be tuned
    smax = floor(N/4.0);    % approx, can be tuned
    svec = smin:smax;    
    
    F2 = cell(1, length(s));
    s_counter = 1;
    for s = svec % loop over sizes
        Ns = floor(N/s); % number of segments
        Fvws2 = zeros(Ns);
        for v = 1:Ns % loop over segments
            ivec = (v-1)*s + (1:s);

            % Cumulative sum
            u = cumsum(X(ivec));

            % Fit to different polynomials
            x = ivec';  % transpose to get column vector
            y = u;      % already column vector

            % Select polynomial
            % % eq (2), ai + bj +c
            A = [x ones(size(x))];
            % % eq (3), ai^2 + bj^2 + c
%             A = [x.*x ones(size(x))];

            % Find coefficients using solve(A,B) == A\B
            coeff = A\y;
            U = A*coeff;

            % Find residual matrix eps_{u,v}
            eps = u - U;

            % Detrended fluctuation function F^2(u,w,s), eq (8)
            Fvws2(v) = sum(sum(eps.*eps))/(s^2);
        end
        F2{s_counter} = Fvws2;
        s_counter = s_counter + 1;
    end
    
    svec = smin:smax;
    
    % Overall detrended fluctuation, eq. (9)
    F2 = cellfun(@mean2, F2);

    x = log10(svec);
    y = log10(sqrt(F2));
    p = polyfit(x,y,1);
    H = p(1);
    fprintf('H = %1.4f\n', H);

    figure;
    plot(x, y);
    hold all;
    x = log10([smin/2 smax*2]);
    y = polyval(p, x);
    plot(x, y);
    legend('log10(sqrt(FF))', sprintf('%1.4fx + %1.4f', p(1), p(2)), 'Location', 'Best');
end
