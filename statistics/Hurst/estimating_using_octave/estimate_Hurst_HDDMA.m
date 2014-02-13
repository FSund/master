function [H, n_vec, sigma_DMA_squared] = estimate_Hurst_HDDMA(f, theta, n_min, n_max, n_step, show_plot)

%
% Higher-Dimensional Detrending Moving Average
% See "Algorithm to estimate the Hurst exponent of high-dimensional fractals", Anna Carbone, 2007, http://arxiv.org/abs/0711.2892
% and "Detrending moving average algorithm for multifractals", Gao-Feng Gu and Wei-Xing Xhou, 2010, http://arxiv.org/abs/1005.0877v2
%
%
% The procedure [H, nvec, sigma_DMA_squared] = estimate_Hurst_HDDMA(f, n_max, theta, show_plot)
% is used to calculate the Hurst exponent of two-dimensional data (for example a heightmap).
%
% Input:
%   f: the two-dimensional data
%   n_min: the lower bound of the segment size n
%   n_max: the upper bound of the segment size n
%   theta: the position parameter of the moving window
%
% Output:
%   H: Hurst exponent
%   nvec: vector with segment sizes
%   sigma_DMA_squared: variance
%
% %% FROM Gu2010 source code (se arxiv.org) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The procedure works as follows:
%   1) For each n, construct the cumulative sum Y in a moving window.
%   2) Calculate the moving average function \widetilde{Y}.
%   3) Determine the residual e by detrending \widetilde{Y} from Y.
%   4) Estimate the root-mean-square function F.
%   5) Calculate the q-th order overall fluctuation function Fq.
%   6) Calculate the multifractal scaling exponent tau(q).
%   7) Calculate the singularity strength function alpha(q) and spectrum f(alpha).
%
% Note:
%   2) The lower bound n_min would better be selected around 10.
%   3) The upper bound n_max would better be selected around 10% of min(size(X)).
%   4) N would better be seleceted in the range [20,40].
%   5) The parameter theta varies in the range [0,1], and we have
%      theta=theta_1=theta_2. Theta = 0 corresponds to backward MFDMA, and
%      theta = 0.5 corresponds to the centered MFDMA, and theta = 1 corresponds to
%      the forward MFDMA. We recommend theta=0.
%   6) In the procedure, we have n=n_1=n_2 for the segment size.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example:
% 1D:
%   [H, nvec, sigma_DMA_squared] = estimate_Hurst_HDDMA(rand(1,1000), 10, 20, 2, 0.0)
% 2D:
%   [H, nvec, sigma_DMA_squared] = estimate_Hurst_HDDMA(rand(100,100), 2, 10, 2, 0.0)
%

if nargin < 6
    show_plot = false;
end

if ( (nargout == 0) && (nargin < 6) )
    show_plot = true;
end

if (nargin > 2 && nargin < 5 || nargin == 0 || (nargin < 2 && nargin > 0))
    error('Too few input arguments')
end

s = size(f);
N = min(s(s > 1));

if (nargin == 2)
    n_min = 2;
    n_max = floor(N/4);
    n_step = 1;
end

dim = sum(size(f) > 1);
if (dim > 2)
    error('Input dimension too high');
end
if (dim < 1)
    error('Input dimension too low');
end
if (n_max <= 2)
    error('n_max should be > 2')
end
if (n_max < n_min)
    error('n_max should be > n_min')
end
if (theta > 1.0)
    error('theta should be in [0.0, 1.0]')
end

n_vec = n_min:n_step:n_max;
sigma_DMA_squared = zeros(1, length(n_vec));

% debug %
n_min
n_max
n_step
theta
show_plot
% debug %

if dim==1
    
    if false
        % old implementation
        for k = 1:length(nvec)
            n = nvec(k);

            % See Gu2010, equation (2) for these ranges
            m_lower = ceil((n-1)*(1-theta));
            m_upper = -floor((n-1)*theta);

            sigma_DMA_squared(k) = 0;
            for i = (m_lower+1):(N+m_upper) % +1 because matlab starts at 1
                imax = i - m_upper;
                imin = i - m_lower;
                moving_average = mean( f(imin:imax) );
                sigma_DMA_squared(k) = sigma_DMA_squared(k) + (f(i) - moving_average)^2;
            end
            sigma_DMA_squared(k) = sigma_DMA_squared(k)/(N-n_max); % See Carbone2007, equation (6)
        end
    end
    
    if true
        % new implementation
        for i_n = 1:length(n_vec) % loop over window sizes
            n = n_vec(i_n);
            m = floor(n*theta);
            
            if theta == 1.0
                % If theta == 1.0, then m == n, so (n-m) == 0, and f(0) is
                % illegal in Matlab. So we add +1 to the start point
                start = n - m + 1;
            else 
                start = n - m;
            end
            for i = start:(N-m)
                % Manual -- faster than sum()/n and mean()
                f_tilde_sum = 0.0;
                for k = (i-n+1+m):(i+m) % reverse order compared to article, to get ascending numbers
                    f_tilde_sum = f_tilde_sum + f(k);
                end
                sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n) + (f(i) - f_tilde_sum/n)^2;

%                 % Auto - Slower than manual
%                 f_tilde = mean(f((i-n+1+m):(i+m)));
%                 f_tilde = sum(f((i-n+1+m):(i+m)))/n;
%                 sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n) + (f(i) - f_tilde)^2;
            end
            sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n)/(N - max(n_vec));
        end
    end
end

if dim==2

    if false
        % Old implementation
        for k = 1:length(n_vec) % loop over window sizes, using square windows (n1 == n2)
            n = n_vec(k);

            % See Gu2010, equation (2) for these ranges
            m_lower = ceil((n-1)*(1-theta));
            m_upper = -floor((n-1)*theta);

            sigma_DMA_squared(k) = 0;
            for i = (m_lower+1):(N+m_upper) % loop over window x-positions
                imax = i - m_upper;
                imin = i - m_lower;
                for j = (m_lower+1):(N+m_upper) % loop over window y-positions
                    % the window is imin:imax, jmin:jmax
                    jmax = j - m_upper;
                    jmin = j - m_lower;
                    moving_average = mean(mean( f(imin:imax, jmin:jmax) )); % average of window

                    % NOTE: There is a difference between mean(mean()) as in moving_average above, and sum(sum())/n^2, but 
                    %       the error is of the size e-16. mean(mean()) is probably faster!
                    % Alternative is: moving_average = sum(sum(f(imin:imax, jmin:jmax)))/(n^2);

                    sigma_DMA_squared(k) = sigma_DMA_squared(k) + (f(i,j) - moving_average)^2;
                end
            end
            sigma_DMA_squared(k) = sigma_DMA_squared(k)*(1/(N-n_max)^dim); % See Carbone2007, equation (9)
        end
    end
    
    if true
        % New implementation
        for i_n = 1:length(n_vec) % loop over window sizes
            n = n_vec(i_n);
            m = floor(n*theta);
            
            if theta == 1.0
                % If theta == 1.0, then m == n, so (n-m) == 0, and f(0) is
                % illegal in Matlab. So we add +1 to the start point
                start = n - m + 1;
            else 
                start = n - m;
            end
            for i = start:(N-m)
                for j = start:(N-m)
                    
%                     % Manual
%                     f_tilde_sum = 0.0;
%                     for k = (i-n+1+m):(i+m) % reverse order compared to article, to get ascending numbers
%                         for l = (j-n+1+m):(j+m)
%                             f_tilde_sum = f_tilde_sum + f(k,l);
%                         end
%                     end
%                     sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n) + (f(i,j) - f_tilde_sum/n^2)^2;

                    % Auto
%                     f_tilde = mean(mean( f((i-n+1+m):(i+m), (j-n+1+m):(j+m)) ));
                    f_tilde = sum(sum( f((i-n+1+m):(i+m), (j-n+1+m):(j+m)) ))/n^2; % fastest
                    sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n) + (f(i,j) - f_tilde)^2;
                end
                sigma_DMA_squared(i_n) = sigma_DMA_squared(i_n)/(N - max(n_vec))^2;
            end
        end
    end
end
        
x = log(dim.*(n_vec.*n_vec));
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