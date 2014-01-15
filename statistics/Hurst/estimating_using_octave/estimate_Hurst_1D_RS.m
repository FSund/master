
% .m

function [H, nvec, rescaled_range] = estimate_Hurst_1D_RS(f, max_n_steps, show_plot)

%
% Estimate the Hurst exponent using rescaled range analysis: http://en.wikipedia.org/wiki/Hurst_exponent
%

if nargin < 3
    show_plot = false;
end

N = length(f);

n_steps = max_n_steps;
% decreasing number of steps until the shortest partial series length is larger than 1
while (2^n_steps >= N)
    n_steps = n_steps-1;
end
nvec = floor(N./(2.^(0:n_steps-1)));
rescaled_range = zeros(1, length(nvec));
for k = 1:length(nvec) % loop over length of partial series
    n = nvec(k);
    for i = 1:n:(N-n+1) % loop over all non-overlapping partial series of length n, in steps of n
        X = f(i:(i+n-1));
        m = mean(X);
        Y = X - m;
        Z = cumsum(Y);
        R = max(Z) - min(Z);
        % S = std(X);
        S = sqrt(sum(Y.^2)/n);  % Uncorrected sample standard deviation
                                % std(X) gives corrected sample standard deviation (same as S above, but replace n with n-1)
        if (S~=0)
            rescaled_range(k) = rescaled_range(k) + R/S;
        end
    end
    rescaled_range(k) = rescaled_range(k)/(N/n); % divide by number of partial series of length n
end

x = log10(nvec);
y = log10(rescaled_range);
fit = polyfit(x, y, 1);
H = fit(1);

if show_plot
    figure;
    plot(x, y);
    hold all;
    plot(x, polyval(fit, x), 'r-');
    legend(['H = ' num2str(H)]);
end
