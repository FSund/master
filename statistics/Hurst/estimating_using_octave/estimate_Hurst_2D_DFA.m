function [H, svec, F2] = estimate_Hurst_2D_DFA(X)

    [M,N] = size(X);
    if (M == 1 || N == 1)
        error('Input needs to be a 2d matrix. Aborting!')
    end
    smin = 6;                       % approx, can be tuned
    smax = floor(min([M,N])/4.0);   % approx, can be tuned
%     svec = smin:smax;
    svec = round(logspace(log10(smin), log10(smax)));
    svec = unique(svec);    % remove duplicates
    
%     % Debug %
%     svec = smax;
%     figure;
%     [m, n] = size(X);
%     [xx, yy] = meshgrid(1:m, 1:n);
%     surface(xx, yy, X, 'FaceAlpha' , 0.5, 'LineStyle', 'none', 'FaceColor', 'r');
%     view(35,10);
%     hold all;
%     default_axis = axis;
%     default_axis(1:4) = [0 m 0 n];
%     % Debug %
    
    F2 = zeros(1, length(svec));
    for is = 1:length(svec) % loop over sizes
        s = svec(is);
        Ms = floor(M/s); % number of segments
        Ns = floor(N/s); % number of segments
        F2sum = 0.0;
        for v = 1:Ms % loop over segments
            for w = 1:Ns % loop over segments
                ivec = (v-1)*s + (1:s);
                jvec = (w-1)*s + (1:s);
                Xuv = X(ivec, jvec);

%                 % Debug %
%                 x = [ivec(1) ivec(1) ivec(end) ivec(end)];
%                 y = [jvec(1) jvec(end) jvec(1) jvec(end)];
%                 z = X(x, y);
%                 line(x, y, z, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'Color', 'r')
%                 % Debug %
                
                % Cumulative sum
                uvw = cumsum(cumsum(X(ivec, jvec),2),1);
                % uvw = cumsum2(Xuv, s);
                % any(abs(uvw - cumsum(cumsum(Xuv,2),1)) > 1e-10)
                
                % Fit to different polynomials
                [xx,yy] = ndgrid(ivec, jvec); % correct!
                x = xx(:);
                y = yy(:);
                z = uvw(:);
                
                % Select polynomial
%                 A = [x      y                       ones(size(x))];
%                 A = [x.*x   y.*y                    ones(size(x))];
%                 A = [x.*y   x       y               ones(size(x))];
                A = [x.*x   y.*y    x       y       ones(size(z))];
%                 A = [x.*x   y.*y    x.*y    x   y   ones(size(z))];
                
                % Find coefficients using solve(A,B) == A\B
                coeff = A\z;
                U = A*coeff;
                
                % Find residual matrix eps_{u,v}
                eps = z - U; % vector
                
                % Debug %
%                 if (v == w)
%                     coeff
%                     [xx, yy] = ndgrid(ivec, jvec);
%                     line(x, y, z, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'Color', 'r')
%                     surface(xx, yy, uvw, 'FaceAlpha' , 0.5, 'FaceColor', 'g', 'LineStyle', 'none');
%                     surface(xx, yy, reshape(U, size(uvw)), 'FaceAlpha' , 0.5, 'FaceColor', 'm', 'LineStyle', 'none');
%                     surface(xx, yy, reshape(eps, size(uvw)), 'FaceAlpha' , 0.5, 'FaceColor', 'y', 'LineStyle', 'none');
%                     surface(xx, yy, zeros(size(uvw)), 'FaceAlpha' , 0.5, 'FaceColor', 'k', 'LineStyle', 'none');
%                 end
                % Debug %

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
%                 Fvws2(v, w) = sum(sum(eps.*eps))/(s^2);
                F2sum = F2sum + sum(eps.*eps)/(s^2); % for vector eps
            end
        end
        % Overall detrended fluctuation, eq. (9)
        F2(is) = F2sum/(Ms*Ns);
    end

%     % Debug %
%     view(35, 10);
% %     axis(default_axis);
%     H = 2.0;
%     return
%     % Debug %

    x = log10(svec);
    y = log10(sqrt(F2));
    p = polyfit(x, y, 1);
    H = p(1);
%     fprintf('H = %1.4f\n', H);

%     figure;
%     plot(x, y, '.');
%     hold all;
%     plot(log10([6 2e2]), log10([1.5e-1 2e0]));
%     plot(log10([6 2e2]), log10([1.5e-1 2e0]*10));
%     axis([5 6e2 1e-1 5e1]);
%     legend('log10(sqrt(FF))', 'article');


%     function uvw = cumsum2(Xuv, s)
%         uvw = zeros(s,s);
%         for i=1:s
%             for j=1:s
%                 cumsum = 0.0;
%                 for k1=1:i
%                     for k2=1:j
%                         cumsum = cumsum + Xuv(k1, k2);
%                     end
%                 end
%                 uvw(i,j) = cumsum;
%             end
%         end
%     end
end
