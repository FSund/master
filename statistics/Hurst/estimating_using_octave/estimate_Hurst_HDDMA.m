
% estimate_Hurst_HDDMA.m

function [H, nvec, sigma_DMA_squared] = estimate_Hurst_HDDMA(f, n_max, theta, show_plot)

%
% Higher-Dimensional Detrending Moving Average
% See "Algorithm to estimate the Hurst exponent of high-dimensional fractals", Anna Carbone, 2007, http://arxiv.org/abs/0711.2892
% and "Detrending moving average algorithm for multifractals", Gao-Feng Gu and Wei-Xing Xhou, 2010, http://arxiv.org/abs/1005.0877v2
%

if nargin < 4
    show_plot = false;
end

dim = sum(size(f) > 1);
if (dim > 2)
    disp('Input dimension too high');
    return;
end
if (dim < 1)
    disp('Input dimension too high');
    return;
end
if (n_max <= 2)
    disp('n_max should be > 2')
    return;
end
if (theta > 1.0)
    disp('theta should be < 1.0')
    return;
end


N = length(f);
if (n_max/N > 0.1)
    disp(['n_max should be << N = ' num2str(N)]);
    return;
end
n_min = 2;
nvec = n_min:n_max;

sigma_DMA_squared = zeros(1,length(nvec));
if dim==1
    %% 1d %%
    dim
    for k = 1:length(nvec)
        n = nvec(k);

        % See Gu2010, equation (2) for these ranges
        m_lower = ceil((n-1)*(1-theta));
        m_upper = -floor((n-1)*theta);

        sigma_DMA_squared(k) = 0;
        for i = (m_lower+1):(N+m_upper) % +1 because matlab starts at 1
            imax = i - m_upper;
            imin = i - m_lower;
            % i
            % imax
            % imin
            moving_average = mean( f(imin:imax) );
            % moving_average
            sigma_DMA_squared(k) = sigma_DMA_squared(k) + (f(i) - moving_average)^2;
        end
        sigma_DMA_squared(k) = sigma_DMA_squared(k)/(N-n_max);
    end

    x = log(nvec.*nvec);
    y = log(sigma_DMA_squared);
    fit = polyfit(x, y, 1);
    H = fit(1);

    if show_plot
        figure;
        plot(x, y)
        hold all;
        plot(x, polyval(fit, x), 'r-')
        title(['1d, theta = ' num2str(theta)])
        legend(['H = ' num2str(H)])
    end
end

if dim==2
    %% 2d %%
    dim
    for k = 1:length(nvec) % loop over window sizes, using square windows (n1 == n2)
        n = nvec(k);
        m = floor(n*theta);
        if (theta == 1.0)
            start = n-m+1;
        else
            start = n-m;
        end
        sigma_DMA_squared(k) = 0;
        for i = start:(N-m) % loop over window positions
            imax = i+m;
            imin = imax - n+1;
            for j = start: (N-m) % loop over window positions
                % the window is imin:imax, jmin:jmax
                jmax = j+m;
                jmin = jmax - n+1;
                moving_average = mean(mean( f(imin:imax, jmin:jmax) )); % average of window

                % NOTE: There is a difference between mean(mean()) as in moving_average above, and sum(sum())/n^2, but 
                %       the error is of the size e-16. mean(mean()) is probably faster!
                % moving_average = sum(sum(f(imin:imax, jmin:jmax)))/(n^2);

                sigma_DMA_squared(k) = sigma_DMA_squared(k) + (f(i,j) - moving_average)^2;
            end
        end
        sigma_DMA_squared(k) = sigma_DMA_squared(k)*(1/(N-n_max)^dim);
    end

    x = log(dim.*(nvec.*nvec));
    y = log(sigma_DMA_squared);
    fit = polyfit(x, y, 1);
    H = fit(1);

    if show_plot
        figure;
        plot(x, y)
        hold all;
        plot(x, polyval(fit, x), 'r-')
        title(['2d, theta = ' num2str(theta)])
        legend(['H = ' num2str(H)])
    end
end
