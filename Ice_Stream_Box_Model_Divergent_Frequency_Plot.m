%% Compute and Plot Frequency Comparison (Grimstead et al., ,2026)
% This is the plot code for the high resolution parameter sweep 
% for the divergent ice stream box model described in the MSc dissertation
%  of Joshua Grimstead and subsequently Grimstead et al., ,2026.  
%
% In this script, the Ice_Stream_Box_Model_Divergent_Parameter_Sweep.m
% output is plotted in three different ways, by comparing the frequencies
% of box 2 and box 3's void ratios. 
%  
%
% 
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5. 
%
%%
load('raw_simulation2_output3_HPC.mat', ...
    'n_points_mat', 'n_points_mat_e3', ...
    'PC_results_e2', 'PC_results_e3', ...
    'peak_times_e2_mat', 'peak_times_e3_mat', ...
    'ThermalDampeningMat', 'FailedPixelsMask', ...
    'T_s2_vals', 'L1_vals');
%%

% 1. Pre-allocate ratio matrix
freq_ratio_mat = NaN(size(n_points_mat));
max_denominator = 32; % Threshold: Higher = more permissive, Lower = only simple locks

% 2. Calculate average frequencies and filter by "rational complexity"
for h = 1:length(T_s2_vals)
    for i = 1:length(L1_vals)
        
        if FailedPixelsMask(h, i) || n_points_mat(h, i) < 10 || n_points_mat_e3(h, i) < 10
            continue; 
        end
        
        t2 = peak_times_e2_mat{h, i};
        t3 = peak_times_e3_mat{h, i};
        
        if ~isempty(t2) && ~isempty(t3) && length(t2) > 2 && length(t3) > 2
            f2 = (length(t2)-1) / (t2(end) - t2(1));
            f3 = (length(t3)-1) / (t3(end) - t3(1));
            
            if f3 > 0
                ratio = f2 / f3;
                
                % Convert to simple fraction and filter
                [p, q] = rat(ratio, 0.005); % 0.005 is the tolerance
                
                % Only assign if the fraction is "simple" (low denominator) 
                % AND within our 0-5 range
                if q <= max_denominator && ratio <= 5
                    freq_ratio_mat(h, i) = ratio;
                end
            end
        end
    end
end

% 3. Visualization
figure('Color', 'w', 'Name', 'Phase Locked Tongues (Complexity Filtered)');
imagesc(T_s2_vals, L1_vals, freq_ratio_mat');
set(gca, 'YDir', 'normal');
caxis([0 3]); 
colorbar;
ylabel('L1'); xlabel('T_{s2}');
title(['Ratio of Frequencies (Filtered: q \leq ' num2str(max_denominator) ')']);
colormap(turbo); 
set(gca, 'Color', [0.8 0.8 0.8]);

%% Complexity and Winding using e2 e3
% 1. Pre-allocate the maps
ComplexityMap = NaN(size(n_points_mat));
WindingMap    = NaN(size(n_points_mat));

% 2. Calculate Rationality based on peak timing ratios
% We use the ratio of (e2 / e3) as the basis for the winding
for h = 1:size(n_points_mat, 1)
    for i = 1:size(n_points_mat, 2)
        
        % Filter: Skip failed/low-count data
        if FailedPixelsMask(h, i) || n_points_mat(h, i) < 4 || n_points_mat_e3(h, i) < 4
            continue; 
        end
        
        t2 = peak_times_e2_mat{h, i};
        t3 = peak_times_e3_mat{h, i};
        
        if ~isempty(t2) && ~isempty(t3) && length(t2) > 2 && length(t3) > 2
            f2 = (length(t2)-1) / (t2(end) - t2(1));
            f3 = (length(t3)-1) / (t3(end) - t3(1));
            
            if f3 > 0
                ratio = f2 / f3;
                W = mod(ratio, 1); % Wrap to [0, 1]
                
                % Use rat() to find the denominator (complexity)
                % We look for simple fractions (p/q) within 1% tolerance
                [p, q] = rat(W, 0.01);
                
                % Complexity Score: 1/q (Simple fractions = High Score)
                ComplexityMap(h, i) = 1 / q;
                WindingMap(h, i)    = W;
            end
        end
    end
end

% 3. Visualization: The "Complexity" Plot (Arnold Tongues)
figure('Color', 'w', 'Name', 'Arnold Tongues: Complexity Score');
ax1 = axes;
h1 = imagesc(T_s2_vals,L1_vals,  ComplexityMap');
set(h1, 'AlphaData', ~isnan(ComplexityMap')); % Punch holes in chaos
set(ax1, 'Color', [0.2 0.2 0.2]); 
colormap(ax1, parula);
axis xy; colorbar;
title('Arnold Tongues (1/Denominator)');
ylabel('L1'); xlabel('T_{s2}');

% 4. Visualization: The Winding Number Plot
figure('Color', 'w', 'Name', 'Arnold Tongues: Winding Number');
ax2 = axes;
h2 = imagesc( T_s2_vals, L1_vals, WindingMap');
set(h2, 'AlphaData', ~isnan(WindingMap'));
set(ax2, 'Color', 'k'); % Black background makes HSV colors pop
colormap(ax2, hsv); 
axis xy; colorbar;
title('Winding Number (Relative Phase)');
ylabel('L1'); xlabel('T_{s2}');

%% 

% Define the output filename
output_file = 'arnold_tongue_complexity_winding_frequencyratio.mat';

% Save the generated maps alongside the original coordinate vectors
save(output_file, ...
     'ComplexityMap', ...    % 1/q simplicity matrix
     'WindingMap', ...    % Logical matrix where LyapMat < 0
     'freq_ratio_mat', ...
     'T_s2_vals', ...        % X-axis grid coordinates
     'L1_vals');             % Y-axis grid coordinates

fprintf('Saved successfully. File ready for local plotting: %s\n', output_file);