function [lambda, lambda_phys, cloud_rel] = Ice_Stream_Box_Model_Divergent_Rosenstein_Grid(PC_results_e2_all, crossing_times_e2_mat, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rosenstein LLEs, Parameter Sweep (Grimstead et al., 2026)
% This is the full parameter sweep which estimates lyapunov exponents
% by running Ice_Stream_Box_Model_Divergent_Rosenstein_Cell.m for the 
% divergent ice stream box model described in the MSc dissertation of 
% Joshua Grimstead and subsequently Grimstead et al., ,2026.
%
%
% Returns three matrices. lambda is the exponent per
% crossing, i.e. the exponent of the Poincare map itself, lambda_phys
% is that exponent per kyr and is NaN wherever crossing times were missing 
% or invalid. cloud_rel is the normalised attractor size.
%
% See Ice_Stream_Box_Model_Divergent_Overlay_Switches.m for the plotting of
% these outputs.
%
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5.
%
%
%%
[num_T, num_L] = size(PC_results_e2_all);

% Preallocated, unlike the generating script's PC_results_e2_all: a
% parfor sliced output keeps whatever was in unassigned elements, so an
% unpreallocated variable can silently carry stale cells from a previous
% run into a new result.
lambda      = NaN(num_T, num_L);
lambda_phys = NaN(num_T, num_L);
cloud_rel   = NaN(num_T, num_L);

fprintf('Rosenstein grid: %d x %d cells\n', num_T, num_L);
tic

dq = parallel.pool.DataQueue;
count = 0;
afterEach(dq, @(~) printProgress());

parfor h = 1:num_T
    row_lam = NaN(1, num_L);
    row_lph = NaN(1, num_L);
    row_cld = NaN(1, num_L);

    for i = 1:num_L
        PC = PC_results_e2_all{h,i};
        ct = crossing_times_e2_mat{h,i};
        if isempty(PC) || isempty(ct)
            continue
        end
        r = Ice_Stream_Box_Model_Divergent_Rosenstein_Cell(PC, ct, varargin{:});
        row_lam(i) = r.lambda;
        row_lph(i) = r.lambda_phys;
        row_cld(i) = r.cloud_rel;
    end

    lambda(h,:)      = row_lam;
    lambda_phys(h,:) = row_lph;
    cloud_rel(h,:)   = row_cld;

    send(dq, h);
end

fprintf('done in %.1f s\n', toc);

    function printProgress()
        count = count + 1;
        fprintf('Processed Row %d of %d (%.1f%%)\n', count, num_T, (count/num_T)*100);
    end
end