%% Plot Frequency Comparison and LLEs w/ Switches (Grimstead et al., 2026)
% This is the plot code for the high resolution parameter sweep and the
%  Rosenstein algorithm outputs with with non-linear switches overlaid, 
% for the divergent ice stream box model described in the MSc dissertation
% of Joshua Grimstead and subsequently Grimstead et al., ,2026.  
%
% In this script, the Ice_Stream_Box_Model_Divergent_Frequency_Plot.m
% output is plotted, as is the 
% Ice_Stream_Box_Model_Divergent_Rosenstein_Grid.m output, with the
% dynamical switches from Ice_Stream_Box_Model_Divergent_Parameter_Sweep.m
% overlaid. 
%  
%
% 
%
% Code written by Joshua Grimstead, revisions made using Claude Opus 5. 
%
%
%%

% 1. Load all required data matrices
load('raw_simulation2_output3_HPC.mat', 'ThresholdCrossingIndiciesMat', 'T_s2_vals', 'L1_vals');
load('arnold_tongue_complexity_winding_frequencyratio.mat', 'ComplexityMap', 'WindingMap', 'freq_ratio_mat');
Lkm = L1_vals / 1e3;   % L1 in km for the axis
load('rosenstein_grid_output3.mat', 'lambda', 'cloud_rel', 'lambda_phys');
%%
% define events
event_names = {
    'Thermal Dampening IceStream 1', ...
    'Thermal Dampening IceStream 2', ...
    'Thermal Dampening IceStream 3', ...
    'Max Till IceStream 1', ...
    'Max Till IceStream 2', ...
    'Max Till IceStream 3', ...
    'Min Void Ratio IceStream 1', ...
    'Min Void Ratio IceStream 2', ...
    'Min Void Ratio IceStream 3'
};

%% Extract boundaries
fprintf('Extracting 9 clean boundary lines...\n');
Boundaries = cell(1, 9);

for k = 1:9
    EventMask = double(ThresholdCrossingIndiciesMat(:, :, k));
    C = contourc(L1_vals, T_s2_vals, EventMask, [0.5 0.5]);
    
    idx = 1;
    segments = {};
    while idx < size(C, 2)
        num_pts = C(2, idx);
        x_coords = C(1, idx+1 : idx+num_pts);
        y_coords = C(2, idx+1 : idx+num_pts);
        segments{end+1} = [x_coords(:), y_coords(:)]; 
        idx = idx + num_pts + 1;
    end
    Boundaries{k} = segments;
end
%% Plot without overlay
% Limit L1 to 120e3
L1_cap  = 120e3;                 
keep_L1 = L1_vals <= L1_cap;
L1_plot = L1_vals(keep_L1);

% Backgrounds, colormaps, titles, and a log-scale flag
plot_data = {
    cloud_rel,      'parula', 'Relative Attractor Diameter',            true;
    ComplexityMap,  'parula', 'Arnold Tongues (Simplicity Score: 1/q)', false;
    WindingMap,     'hsv',    'Phase Locking Diagram (Winding Map)',    false;
    freq_ratio_mat, 'jet',    'Frequency Ratio Map',                    false
};

for f = 1:4
    grid_data  = plot_data{f, 1};
    cmap       = plot_data{f, 2};
    plot_title = plot_data{f, 3};
    use_log    = plot_data{f, 4};

    grid_data = grid_data(:, keep_L1);        % apply the L1 cap

    if use_log % use a log colour scale
        grid_data(grid_data <= 0) = NaN;      % log needs strictly positive
    end

    figure('Color', 'w', 'Position', [100 + (f*50), 100, 950, 450]);
    hold on;

    h_img = imagesc(T_s2_vals, L1_plot, grid_data');
    set(gca, 'YDir', 'normal', 'Color', [0.8 0.8 0.8]);
    set(h_img, 'AlphaData', ~isnan(grid_data'));
    colormap(gca, cmap);

    if use_log
        set(gca, 'ColorScale', 'log');
        v = grid_data(isfinite(grid_data));
        if ~isempty(v), clim([min(v) max(v)]); end
    end
    colorbar;

    axis xy tight;
    ylim([min(L1_plot) max(L1_plot)]);   
    xlabel('T_{s2} Values');
    ylabel('L_1 Values');
    title(plot_title);
    legend('Location', 'eastoutside', 'Box', 'off');
    grid on;
end

%% Plot with overlay
% --- L1 cap -----------------------------------------------------------
L1_cap  = 120e3;                
keep_L1 = L1_vals <= L1_cap;
L1_plot = L1_vals(keep_L1);

% Backgrounds, colormaps, titles, and a log-scale flag
plot_data = {
    cloud_rel,      'parula', 'Relative Attractor Diameter',            true;
    ComplexityMap,  'parula', 'Arnold Tongues (Simplicity Score: 1/q)', false;
    WindingMap,     'hsv',    'Phase Locking Diagram (Winding Map)',    false;
    freq_ratio_mat, 'jet',    'Frequency Ratio Map',                    false
};

for f = 1:4
    grid_data  = plot_data{f, 1};
    cmap       = plot_data{f, 2};
    plot_title = plot_data{f, 3};
    use_log    = plot_data{f, 4};

    grid_data = grid_data(:, keep_L1);        % apply the L1 cap

    if use_log
        grid_data(grid_data <= 0) = NaN;      
    end

    figure('Color', 'w', 'Position', [100 + (f*50), 100, 950, 450]);
    hold on;

    h_img = imagesc(T_s2_vals, L1_plot, grid_data');
    set(gca, 'YDir', 'normal', 'Color', [0.8 0.8 0.8]);
    set(h_img, 'AlphaData', ~isnan(grid_data'));
    colormap(gca, cmap);

    if use_log
        set(gca, 'ColorScale', 'log');
        v = grid_data(isfinite(grid_data));
        if ~isempty(v), clim([min(v) max(v)]); end
    end
    colorbar;

    % Overlay only the active boundary lines
    for i = 1:num_active
        original_k = active_indices(i);    % Events
        current_color = active_colors(i, :); % corresponding colours
        segments = Boundaries{original_k};

        for seg_idx = 1:length(segments)
            line_coords = segments{seg_idx};

            x_plot = line_coords(:, 2);
            y_plot = line_coords(:, 1);

            if seg_idx == 1
                plot(x_plot, y_plot, 'Color', current_color, 'LineWidth', 2.5, ...
                     'DisplayName', event_names{original_k});
            else
                plot(x_plot, y_plot, 'Color', current_color, 'LineWidth', 2.5, ...
                     'HandleVisibility', 'off');
            end
        end
    end
    axis xy tight;
    ylim([min(L1_plot) max(L1_plot)]);  
    xlabel('T_{s2} Values');
    ylabel('L_1 Values');
    title(plot_title);
    legend('Location', 'eastoutside', 'Box', 'off');
    grid on;
end

%% Rosenstien Plots
markPeriod1  = true;
period1_val  = 0;

lam  = lambda;
lamp = lambda_phys;
if markPeriod1
    lam(regime == -1)  = period1_val;
    lamp(regime == -1) = period1_val;
end

% --- lambda (per crossing) ------------------------------------------
figure('Color','w','Name','lambda');
imagesc(T_s2_vals, Lkm, lam', 'AlphaData', ~isnan(lam'));
set(gca,'YDir','normal','Color',[.85 .85 .85]);
colorbar; colormap(turbo);
xlabel('T_{s2} (-°C)'); ylabel('L1 (km)');
%title('\lambda per crossing');

% --- lambda_phys (per kyr) ------------------------------------------
figure('Color','w','Name','lambda phys');
imagesc(T_s2_vals, Lkm, lamp', 'AlphaData', ~isnan(lamp'));
set(gca,'YDir','normal','Color',[.85 .85 .85]);
colorbar; colormap(turbo);
xlabel('T_{s2} (-°C)'); ylabel('L1 (km)');
%title('\lambda per kyr');

%% --- e-folding time (predictability time) ---------------------------

tef = 1 ./ lambda_phys;
diverging = lambda_phys > 0;

v  = tef(diverging & isfinite(tef));
lo = prctile(v, 2);
hi = prctile(v, 98);

tef_disp = tef;
tef_disp(~diverging) = hi * 10;      % not diverging -> infinite e-folding

% Apply the period 1 mask
if markPeriod1
    tef_disp(regime == -1) = hi * 10;
end

tef_disp(tef_disp < lo) = lo;

figure('Color','w','Name','e-folding time');
imagesc(T_s2_vals, Lkm, tef_disp', 'AlphaData', ~isnan(tef_disp'));
set(gca,'YDir','normal','Color',[.85 .85 .85],'ColorScale','log');
clim([lo hi]);
colormap(flipud(hot));
cb = colorbar; ylabel(cb, 'kyr');
cb.Ticks = [lo, sqrt(lo*hi), hi];
cb.TickLabels = {sprintf('%.3g kyr',lo), sprintf('%.3g kyr',sqrt(lo*hi)), ...
                 sprintf('\\geq%.3g kyr',hi)};
xlabel('T_{s2} (-°C)'); ylabel('L1 (km)');
%title('e-folding divergence time (kyr)');
























