%% Rosenstein LLEs, Run Rosenstein Algorithm (Grimstead et al., 2026)
% This runs Ice_Stream_Box_Model_Divergent_Rosenstein_Grid.m over the saved
% parameter sweep output and writes the result to disk, for the divergent
% ice stream box model described in the MSc dissertation of Joshua
% Grimstead and subsequently Grimstead et al., 2026.
%
% Dims is set here rather than left to the default so that the columns
% used are recorded alongside the run.
%
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5.
%
%
%%
dims_used = [1 5 6 7 8 9 10];

if isempty(gcp('nocreate'))
    parpool(48);
end

load('raw_simulation2_output3_HPC.mat', ...
     'PC_results_e2_all', 'crossing_times_e2_mat');

[lambda, lambda_phys, cloud_rel] = ...
    Ice_Stream_Box_Model_Divergent_Rosenstein_Grid( ...
        PC_results_e2_all, crossing_times_e2_mat, 'Dims', dims_used);

save('rosenstein_grid_output3.mat', ...
     'lambda', 'lambda_phys', 'cloud_rel');