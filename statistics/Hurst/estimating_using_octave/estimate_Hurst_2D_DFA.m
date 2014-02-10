function estimate_Hurst_2D_DFA(X)

    [M,N] = size(X);
    smin = 6;                       % approx, can be tuned
    smax = floor(min([M,N])/4.0);   % approx, can be tuned
    F2 = cell(1, length(smin:smax));
    svec = smin:smax;

    s_counter = 1;
    for s = svec % loop over sizes
        Ms = floor(M/s); % number of segments
        Ns = floor(N/s); % number of segments
        Fvws2 = zeros(Ms, Ns);
        for v = 1:Ms % loop over segments
            for w = 1:Ns % loop over segments
                ivec = (v-1)*s + (1:s);
                jvec = (w-1)*s + (1:s);
%                 Xuv = X(ivec, jvec); % not needed
                
                % Cumulative sum
                u = cumsum(cumsum(X(ivec, jvec),2),1); % 1 in the end is default
                
                % Fit to different polynomials              
                [xx,yy] = meshgrid(ivec, jvec);
                x = xx(:); % Reshape/create 1d vectors
                y = yy(:);
                z = u(:);
                
                % Select polynomial
                % % eq (2), ai + bj +c
%                 A = [x y ones(size(x))]
                % % eq (3), ai^2 + bj^2 + c
                A = [x.*x y.*y ones(size(x))];
                
                % Find coefficients using solve(A,B) == A\B
                coeff = A\z;
                U = reshape(A*coeff, size(u));
                
                % Find residual matrix eps_{u,v}
                eps = u - U;

                % % Debug stuff % %
%                 figure(3)
%                 hold on;
%                 line(x, y, z, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'Color', 'r')
%                 surface(xx, yy, U, 'FaceAlpha' , 0.5)
%                 surface(xx, yy, eps, 'FaceAlpha', 0.5)
% %                 surface(xx, yy, X(ivec, jvec), 'FaceAlpha', 0.5)
%                 shading flat;
%                 grid on;           
%                 return
                % % Debug stuff % %

                % Detrended fluctuation function F^2(u,w,s), eq (8)
%                 Fvws(v, w) = sqrt( sum(sum(eps.*eps)) / (s^2) );
%                 Fvws2(v, w) = mean2(eps.*eps);
                Fvws2(v, w) = sum(sum(eps.*eps))/(s^2);
            end
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
    fprintf('H = %1.2f\n', p(1));

    figure;
    plot(x, y);
    hold all;
    plot(log10([6 2e2]), log10([1.5e-1 2e0]));
%     plot(log10([6 2e2]), log10([1.5e-1 2e0]*10));
%     axis([5 6e2 1e-1 5e1]);
    legend('log10(sqrt(FF))', 'article');
end
