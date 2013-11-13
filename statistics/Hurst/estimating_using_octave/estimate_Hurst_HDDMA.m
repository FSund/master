
% GFGU_MFDMA_2D.m

function [H, nvec, sigma_DMA_squared] = estimate_Hurst_HDDMA(f, n_max, theta)

%
% Higher-Dimensional Detrending Moving Average
%


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
    disp('n_max should be << N');
    return;
end
n_min = 2;
nvec = n_min:n_max;

sigma_DMA_squared = zeros(1,length(nvec));
if dim==1
    %% 1d %%
    for n = nvec
        m = floor(n*theta); % integer part of positive number
        if (theta == 1.0)
            start = n-m+1;
        else
            start = n-m;
        end
        sigma_DMA_squared(n-n_min+1) = 0;
        for i = start:(N-m) % +1 because matlab starts at 1
            imax = i+m;
            imin = imax - n+1;
            moving_average = mean( f(imin:imax) );
            sigma_DMA_squared(n-n_min+1) = sigma_DMA_squared(n-n_min+1) + (f(i) - moving_average)^2;
        end
        sigma_DMA_squared(n-n_min+1) = sigma_DMA_squared(n-n_min+1)*(1/(N-n_max)^dim); % n-n_min+1 to start at first index in vector
    end

    x = log(dim.*(nvec.*nvec));
    y = log(sigma_DMA_squared);

    figure;
    plot(x, y)
    hold all;
    fit = polyfit(x, y, 1);
    H = fit(1);
    plot(x, polyval(fit, x), 'r-')
    title(['1d, theta = ' num2str(theta)])
    legend(['H = ' num2str(H)])
end

if dim==2
    %% 2d %%
    for n = nvec % loop over window sizes, using square windows (n1 == n2)
        m = floor(n*theta);
        if (theta == 1.0)
            start = n-m+1;
        else
            start = n-m;
        end
        sigma_DMA_squared(n-n_min+1) = 0;
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

                sigma_DMA_squared(n-n_min+1) = sigma_DMA_squared(n-n_min+1) + (f(i,j) - moving_average)^2;
            end
        end
        sigma_DMA_squared(n-n_min+1) = sigma_DMA_squared(n-n_min+1)*(1/(N-n_max)^dim);
    end

    x = log(dim.*(nvec.*nvec));
    y = log(sigma_DMA_squared);

    figure;
    plot(x, y)
    hold all;
    fit = polyfit(x, y, 1);
    H = fit(1);
    plot(x, polyval(fit, x), 'r-')
    title(['2d, theta = ' num2str(theta)])
    legend(['H = ' num2str(H)])
end
