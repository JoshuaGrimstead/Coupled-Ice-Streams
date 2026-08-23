function res = Ice_Stream_Box_Model_Divergent_Rosenstein_Cell(PC_all, cross_times, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rosenstein LLEs, Cell (Grimstead et al., 2026)
% This is the adapted Rosenstein algorithm used to estimate the Largest
% Lyapunov Exponent at a single parameter value, used for the
% divergent ice stream box model described in the MSc dissertation of
% Joshua Grimstead and subsequently Grimstead et al., ,2026.
%
% In this function, a single 12D Poincare sections produced by
% Ice_Stream_Box_Model_Divergent_Parameter_Sweep.m is used to compute a
% lyapunov exponent. Nearest-neighbour pairs are tracked forward through
% the section sequence and the slope of their mean log-separation L(k) is
% fitted to give the exponent. The fit window is set by the size of the 
% normalised attractor to control relative saturation. he fit begins at 
% KMin = 1 to exclude the first step of L(k), as it reflects initial
% nearest neighbour distance and not seperation.
% 
%
% Returns lambda (per crossing), lambda_phys (per kyr, via the mean
% return time), and attractor size - cloud_rel. 
%
% See Ice_Stream_Box_Model_Divergent_Rosenstien_Grid.m for the implemtation
% across the full parameter spread.
%
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5.
%
%
%% Define inputs
ip = inputParser;
ip.addParameter('Dims', [1 5 6 7 8 9 10]);
ip.addParameter('TimeScale', 3.15e10);
ip.addParameter('DegenerateTol', 1e-8);
ip.addParameter('KmaxFrac', 0.25);
ip.addParameter('SatFrac', 0.5);
ip.addParameter('MinFitPts', 5);

ip.addParameter('KMin', 1);
ip.parse(varargin{:});
o = ip.Results;
res = struct('lambda', NaN, 'lambda_phys', NaN, 'cloud_rel', NaN, 'note', '');

%% Validate Poincare section
if isempty(PC_all) || size(PC_all,1) < 15
    res.note = 'fewer than 15 crossings';
    
    return
end
N = size(PC_all,1);


t = cross_times(:) / o.TimeScale;
if numel(t) ~= N
    warning('rosensteinCell:mismatch', ...
        'PC_all has %d rows but %d crossing times; lambda_phys will be NaN.', ...
        N, numel(t));
    t = [];
end

if numel(t) >= 2; mrt = mean(diff(t)); else; mrt = NaN; end


%% 1. Attractor size normalised
% Per column: spread relative to that column's own magnitude. This is
% the only way to compare an ice thickness (~1e2) against a void ratio
% (~1e-1) against whatever sits at ~1e-5.
X = PC_all(:, o.Dims);
sd    = std(X, 0, 1);
scale = max(abs(X), [], 1);          % robust when a column is ~zero-mean
scale(scale == 0) = 1;
rel   = sd ./ scale; % relative attractor size 

keep = rel > o.DegenerateTol;


if ~any(keep)
   
    res.cloud_rel = max(rel);
    res.note = 'every retained column is constant to solver precision';
    return
end

res.cloud_rel = max(rel(keep));


%% 2. Normalise and measure the attractor point cloud
Xk = X(:, keep);
Z  = (Xk - mean(Xk,1)) ./ std(Xk, 0, 1);

D = pdist2(Z, Z);
D(1:N+1:end) = Inf;

theiler = getTheiler(Z);
for r = 1:N
    D(r, max(1,r-theiler):min(N,r+theiler)) = Inf;
end

finite_d = D(isfinite(D));
if numel(finite_d) < 30
    res.note = 'too few admissible pairs after Theiler exclusion';
    return
end
% 90th percentile rather than the max: one outlying crossing shouldn't
% set the saturation ceiling for the whole cell.
cloud_size = prctile(finite_d, 90);


%% 3. Track trajectory pairs
[dmin, jn] = min(D, [], 2);
valid = isfinite(dmin);
if sum(valid) < 10
    res.note = 'too few points with a valid neighbour';
    return
end

kmax = max(5, min(round(o.KmaxFrac * N), N - 1));
logd = NaN(N, kmax + 1);
for k = 0:kmax
    ii = find(valid & (1:N)' + k <= N & jn + k <= N);
    if isempty(ii); continue; end
    dk = vecnorm(Z(ii + k, :) - Z(jn(ii) + k, :), 2, 2);
    dk(dk == 0) = NaN;
    logd(ii, k+1) = log(dk);
end

L = mean(logd, 1, 'omitnan')';   

%% 4. Fit window
sat_level = log(o.SatFrac * cloud_size);

%% 5. Fit LE

k_fit_end = firstAtOrAbove(L, sat_level, o.KMin + o.MinFitPts);
res.lambda = fitSlope(L, k_fit_end, o.MinFitPts, o.KMin);
res.lambda_phys = res.lambda / mrt;

end 


%% ===================================================================
function k = firstAtOrAbove(Lc, level, min_pts)
% First index where the curve reaches the saturation level, floored at
% min_pts and capped at the curve length. Callers pass KMin + MinFitPts
% as the floor, so that dropping the first KMin points still leaves
% MinFitPts to fit. Without that, a cell whose L(k) saturates
% immediately -- a low-period orbit -- ends up with too few points after
% the KMin cut and returns NaN instead of the near-zero slope it should.
    k = find(Lc >= level, 1, 'first');
    if isempty(k)
        k = numel(Lc);          % never approaches the diameter: use it all
    end
    k = min(max(k, min_pts), numel(Lc));
end


function [slope, p, kfit, ok] = fitSlope(Lc, k_end, min_pts, k_start)
% Least-squares slope of Lc over k = k_start .. k_end-1, skipping NaNs.
    slope = NaN; p = []; kfit = []; ok = false;
    if nargin < 4 || isempty(k_start); k_start = 0; end
    i0 = k_start + 1;
    if i0 > k_end; return; end
    kk = (k_start:k_end-1)';
    LL = Lc(i0:k_end);
    good = isfinite(LL);
    if sum(good) < min_pts
        return
    end
    kfit = kk(good);
    p = polyfit(kfit, LL(good), 1);
    slope = p(1);
    ok = true;
end


function w = getTheiler(Z)
% Decorrelation lag of the section sequence. Points adjacent in index
% are adjacent in time and would be counted as "close" for that reason
% alone rather than because the orbit returned near itself.
    N = size(Z,1);
    s = vecnorm(Z, 2, 2); s = s - mean(s);
    max_lag = min(N - 2, round(N/4));
    den = sum(s.^2);
    if max_lag < 2 || den == 0
        w = 1; return
    end
    ac = arrayfun(@(L) sum(s(1:end-L).*s(1+L:end))/den, (1:max_lag)');
    idx = find(ac < 1/exp(1), 1, 'first');
    if isempty(idx); idx = max_lag; end
    w = min(max(idx, 1), floor(N/10));
end
