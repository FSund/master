function estimate_Hurst_2D_DFA(X)

%     function ret = cumsum2(X)
%         [m,n] = size(X);
%         ret = zeros(m,n);
%         for i=1:m
%             for j=1:n
%                 sum = 0.0;
%                 for k=1:i
%                     for l=1:j
%                         sum = sum + X(k,l);
%                     end
%                 end
%                 ret(i,j) = sum;
%             end
%         end
%         return
%     end
% 
%     a = rand(10,10);
%     b = cumsum(cumsum(a,2),1)
%     c = cumsum2(a)
%     d = abs(b-c) > 1e-10
% 
%     return

    [M,N] = size(X);
    smin = 6;                       % approx, can be tuned
    smax = floor(min([M,N])/4.0);   % approx, can be tuned
    F2 = cell(1, length(smin:smax));

%     svec = smax/4;
    svec = smin:smax;

    s_counter = 1;
    for s = svec % loop over sizes
        Ms = floor(M/s); % number of segments
        Ns = floor(N/s); % number of segments
        Fvws2 = zeros(Ms, Ns);
        for v = 1:Ms % loop over segments
            for w = 1:Ns % loop over segments
%                 ivec = ((v-1)*s + 1):(v*s);
                ivec = (v-1)*s + (1:s);
%                 jvec = ((w-1)*s + 1):(w*s);
                jvec = (w-1)*s + (1:s);
%                 Xuv = X(ivec, jvec);
                
                % Auto cumulative sum
                u = cumsum(cumsum(X(ivec, jvec),2),1); % 1 in the end is default
%                 u = cumsum(cumsum(X(ivec, jvec),1),2); % equivalent to above
%                 u = cumsum2(X(ivec, jvec));
                
%                 % Manual cumulative sum, eq. (1)
%                 % ! Tested -- Equivalent to the above !
%                 u = zeros(s,s);
%                 x = X(ivec, jvec);
%                 for i=1:s
%                     for j=1:s
% %                         sum_ = 0.0;
% %                         for ii=1:i
% %                             for jj=1:j
% %                                 sum_ = sum_ + x(ii,jj);
% %                             end
% %                         end
% %                         u(i,j) = sum_;
%                         u(i,j) = sum(sum(x(1:i, 1:j)));
%                     end
%                 end
                
                % Fit to different polynomials {\tilde u}
                
                % polynomial: z = a1*x + a2*x^2 + a3*x^2*y + a4*y^2 etc.
                % the coefficients a1,a2 etc. are given by vector A, defined as
                % A=[x x.^2 x.^2.*y y.^2 etc.]\z,
                % so four our problem we get for polynomial in eq (2)
                % u(i,j) = ai + bj + c
                % --> z = ax + by + c
                
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
%                 mean(mean(eps))

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
%         Fvws2
        s_counter = s_counter + 1;
    end

    svec = smin:smax;
    
    % Overall detrended fluctuation, eq. (9)
    % Auto
    F2_auto = cellfun(@mean2, F2);

    % Manual
%     % ! Tested -- F2_mean equal to auto above, within 1e-7 !
%     F_mean = zeros(size(svec));
%     F2_mean = zeros(size(svec));
%     ii = 1;
%     for s = svec
%         Ms = floor(M/s);
%         Ns = floor(N/s);
%         
%         sum_ = 0.0;
%         for v = 1:Ms
%             for w = 1:Ns
%                 sum_ = sum_ + F2{ii}(v,w);
%             end
%         end
%         F_mean(ii) = sqrt(sum_/(Ms*Ns));
%         F2_mean(ii) = sum_/(Ms*Ns);
%         
%         ii = ii+1;
%     end

    FF = F2_auto;

    x = log10(svec);
    y = (1/5.0)*log10(FF);
%     y = log10(sqrt(FF));
    p = polyfit(x,y,1);
    p

    figure;
    plot(log10(svec), log10(sqrt(FF)));
    hold all;
%     loglog(log10(svec), log10(sqrt(sqrt(FF))));
    plot(x, y);
    plot(log10([6 2e2]), log10([1.5e-1 2e0]));
%     plot(log10([6 2e2]), log10([1.5e-1 2e0]*10));
%     axis([5 6e2 1e-1 5e1]);
    title('log10(x) log10(y)')
    legend('log10(sqrt(FF))', 'custom', 'article');
    
    
    
    
    
    
%     figure;
%     loglog(svec, F);
%     axis([5 6e2 1e-1 5e1])
%     x = log10(svec);
%     y = log10(F);
%     figure; plot(x,y);
%     p = polyfit(x,y,1)
%     p = polyfit(y,x,1)
    
    
%     figure;
%     loglog(svec, sqrt(FF));
%     hold all;
%     loglog(svec, sqrt(sqrt(FF)));
%     plot([6 2e2], [1.5e-1 2e0]);
%     plot([6 2e2], [1.5e-1 2e0]*10);
% %     axis([5 6e2 1e-1 5e1]);
%     title('loglog')
% %     legend('sqrt(FF)', 'sqrt(sqrt(FF))', 'article 1', 'article 2');
    
%     figure;
%     plot(svec,sqrt(FF));
%     hold all;
%     plot(svec,sqrt(sqrt(FF)));
%     plot([6 2e2], [1.5e-1 2e0]);
%     plot([6 2e2], [1.5e-1 2e0]*10);
% %     axis([0 6e2 0 5e1]);
%     title('normal')
%     legend('sqrt(FF)', 'sqrt(sqrt(FF))', 'article 1', 'article 2');

%     figure;
%     loglog(svec, F);
%     axis([5 6e2 1e-1 5e1])
%     x = log10(svec);
%     y = 2*log10(FF);
%     figure; plot(x,y);
%     p = polyfit(x,y,1)
%     p = polyfit(y,x,1)
%     axis([5 6e2 1e-1 5e1])
    
    % ret = [ix' ix2' iy' iy2'];
end
