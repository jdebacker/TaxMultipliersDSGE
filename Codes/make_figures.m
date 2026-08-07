% make_figures.m
% =========================================================================
% Generates figures for the paper from posterior Bayesian IRFs.
%
% Requires that Dynare has been run with the bayesian_irf option so that
% oo_.PosteriorIRF.dsge is populated (Median, HPDinf, HPDsup).
%
% For the two multipliers-over-time figures, the estimation var_list must
% include the own-revenue series. Update the .mod file to:
%   estimation(...) y c i h pi R g tr div taxrev
%                   rev_c rev_w rev_k rev_d rev_i spend_ic;
% and rerun Dynare (with bayesian_irf) before running this script.
%
% Files produced (all PNG, 300 DPI, saved to figures/ subfolder):
%   irf_eps_g.png                    -- Response to government spending shock
%   irf_eps_w.png                    -- Response to labor income tax shock
%   irf_eps_k.png                    -- Response to capital income tax shock
%   irf_eps_c.png                    -- Response to consumption tax shock
%   irf_eps_ic.png                   -- Response to investment tax credit shock
%   irf_eps_etau.png                 -- Response to expensing rate shock
%   multipliers_over_time_base.png   -- Method 1: -PV(DY) / PV(D taxrev)
%   multipliers_over_time_rev.png    -- Method 2: -PV(DY) / D(own_rev)_year1
%
% Usage:
%   dynare tax_dsge_est noclearall
%   make_figures
% =========================================================================

if ~exist('M_','var') || ~exist('oo_','var')
    error('Run Dynare first: dynare tax_dsge_est noclearall');
end

if ~isfield(oo_, 'PosteriorIRF') || ~isfield(oo_.PosteriorIRF, 'dsge')
    error(['oo_.PosteriorIRF.dsge not found. Rerun Dynare with the ' ...
           'bayesian_irf option enabled in the estimation() block.']);
end

%% ========================================================================
%  SETUP
%% ========================================================================
% Figure output folder (relative to current working directory)
fig_dir = 'figures';
if ~exist(fig_dir, 'dir')
    [status, msg] = mkdir(fig_dir);
    if status == 0
        error('Could not create folder %s: %s', fig_dir, msg);
    end
    fprintf('Created folder: %s\n', fig_dir);
end

% Shocks to plot (main paper figures)
shocks = {'eps_g', 'eps_w', 'eps_k', 'eps_c', 'eps_ic', 'eps_etau', 'eps_td'};
shock_labels = {'Government Spending', 'Labor Income Tax', ...
                'Capital Income Tax', 'Consumption Tax', ...
                'Investment Tax Credit', 'Expensing Rate', 'Dividend Tax'};

% Response variables per IRF panel (matches Zubairy 2014 Figs 2-3)
resp_vars   = {'y',      'c',           'i',          'h',     'pi',        'R'};
resp_labels = {'Output', 'Consumption', 'Investment', 'Hours', 'Inflation', 'Interest Rate'};

% Own-revenue variable per shock (for Method 2)
own_rev_map = struct( ...
    'eps_g',    'g', ...
    'eps_w',    'rev_w', ...
    'eps_k',    'rev_k', ...
    'eps_c',    'rev_c', ...
    'eps_ic',   'spend_ic', ...
    'eps_etau', 'rev_k', ...
    'eps_td',   'rev_d');

% Figure appearance
band_color = [0.75 0.85 0.95];   % light blue for HPD band
line_color = [0.00 0.30 0.65];   % darker blue for median
zero_color = [0.40 0.40 0.40];   % gray for zero line

Median = oo_.PosteriorIRF.dsge.Median;
HPDinf = oo_.PosteriorIRF.dsge.HPDinf;
HPDsup = oo_.PosteriorIRF.dsge.HPDsup;

%% ========================================================================
%  IRF PANELS (one figure per shock, 6 subplots each)
%% ========================================================================
for s = 1:length(shocks)
    shock = shocks{s};

    % Check all required response fields exist
    missing_fields = {};
    for v = 1:length(resp_vars)
        fname = [resp_vars{v} '_' shock];
        if ~isfield(Median, fname)
            missing_fields{end+1} = fname; %#ok<SAGROW>
        end
    end
    if ~isempty(missing_fields)
        warning('Skipping %s: missing IRF fields (e.g., %s). Was this variable in the estimation var_list?', ...
                shock, missing_fields{1});
        continue;
    end

    fig = figure('Position', [100 100 900 600], 'Visible', 'off', ...
                 'Color', 'white', 'PaperPositionMode', 'auto');

    for v = 1:length(resp_vars)
        fname = [resp_vars{v} '_' shock];
        med   = Median.(fname)(:);
        lo    = HPDinf.(fname)(:);
        hi    = HPDsup.(fname)(:);

        T = length(med);
        t = (1:T)';

        subplot(2, 3, v);
        fill([t; flipud(t)], [lo; flipud(hi)], band_color, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.6);
        hold on;
        plot([1 T], [0 0], '--', 'Color', zero_color, 'LineWidth', 0.5);
        plot(t, med, '-', 'Color', line_color, 'LineWidth', 1.75);
        hold off;

        title(resp_labels{v}, 'FontSize', 11, 'FontWeight', 'bold');
        xlabel('Quarters', 'FontSize', 9);
        xlim([1 T]);
        grid on;
        set(gca, 'FontSize', 9, 'Box', 'on');
    end

    try
        sgtitle(sprintf('Response to a %s Shock', shock_labels{s}), ...
                'FontSize', 13, 'FontWeight', 'bold');
    catch
        annotation('textbox', [0 0.94 1 0.05], 'String', ...
                   sprintf('Response to a %s Shock', shock_labels{s}), ...
                   'HorizontalAlignment', 'center', 'FontSize', 13, ...
                   'FontWeight', 'bold', 'EdgeColor', 'none');
    end

    filename = fullfile(fig_dir, sprintf('irf_%s.png', shock));
    print(fig, filename, '-dpng', '-r300');
    close(fig);
    fprintf('Saved %s\n', filename);
end

%% ========================================================================
%  DISCOUNT FACTOR AND HORIZON SETUP (used by both multiplier figures)
%% ========================================================================
vnames = cellstr(M_.endo_names);
vi_R = find(strcmp(vnames, 'R'));
R_ss = oo_.steady_state(vi_R);

if ~isfield(Median, 'y_eps_g')
    error('Median IRF for y_eps_g missing; cannot build multiplier figures.');
end
T = length(Median.y_eps_g);
disc = R_ss .^ (-(0:T-1))';
Y1_QTRS = 4;

% Colors for the 7 shocks (consistent across both multiplier figures)
line_colors = [
    0.00 0.45 0.74;   % govt spending
    0.85 0.33 0.10;   % labor tax
    0.93 0.69 0.13;   % capital tax
    0.49 0.18 0.56;   % consumption tax
    0.47 0.67 0.19;   % ITC
    0.30 0.75 0.93;   % expensing
    0.64 0.08 0.18    % dividend tax
];

%% ========================================================================
%  MULTIPLIERS OVER TIME -- METHOD 1 (base): -PV(DY) / PV(D taxrev)
%% ========================================================================
fig = figure('Position', [100 100 900 550], 'Visible', 'off', ...
             'Color', 'white', 'PaperPositionMode', 'auto');
hold on;

t = (1:T)';
legend_handles = [];
legend_names   = {};

for s = 1:length(shocks)
    shock = shocks{s};
    y_fld = ['y_' shock];
    if ~isfield(Median, y_fld); continue; end
    y_irf = Median.(y_fld)(:);

    if strcmp(shock, 'eps_g')
        den_fld = 'g_eps_g';
        if ~isfield(Median, den_fld); continue; end
        den_irf = Median.(den_fld)(:);
        mult = cumsum(disc .* y_irf) ./ cumsum(disc .* den_irf);
    else
        den_fld = ['taxrev_' shock];
        if ~isfield(Median, den_fld); continue; end
        den_irf = Median.(den_fld)(:);
        mult = -cumsum(disc .* y_irf) ./ cumsum(disc .* den_irf);
    end

    % Guard against near-zero denominator in first few periods
    cum_den = cumsum(disc .* den_irf);
    mult(abs(cum_den) < 1e-12) = NaN;

    h = plot(t, mult, '-', 'Color', line_colors(s,:), 'LineWidth', 2);
    legend_handles(end+1) = h; %#ok<SAGROW>
    legend_names{end+1}   = shock_labels{s}; %#ok<SAGROW>
end

plot([1 T], [0 0], '--', 'Color', zero_color, 'LineWidth', 0.5);
hold off;

xlabel('Horizon (quarters)', 'FontSize', 11);
ylabel('Present-Value Multiplier', 'FontSize', 11);
title('Fiscal Multipliers by Horizon (Total Revenue Denominator)', ...
      'FontSize', 13, 'FontWeight', 'bold');
legend(legend_handles, legend_names, 'Location', 'best', 'FontSize', 10);
grid on;
xlim([1 T]);
set(gca, 'FontSize', 10, 'Box', 'on');

filename = fullfile(fig_dir, 'multipliers_over_time_base.png');
print(fig, filename, '-dpng', '-r300');
close(fig);
fprintf('Saved %s\n', filename);

%% ========================================================================
%  MULTIPLIERS OVER TIME -- METHOD 2 (rev): -PV(DY) / D(own_rev)_year1
%% ========================================================================
% Check that all own-revenue variables exist in the Bayesian IRF output
missing_own = {};
for s = 1:length(shocks)
    own_var = own_rev_map.(shocks{s});
    fname = [own_var '_' shocks{s}];
    if ~isfield(Median, fname)
        missing_own{end+1} = fname; %#ok<SAGROW>
    end
end

if ~isempty(missing_own)
    warning(['Skipping Method 2 multiplier figure -- missing IRF fields: %s\n' ...
             '  To fix: add these variables to the estimation var_list in the ' ...
             '.mod file:\n' ...
             '    rev_c rev_w rev_k rev_d rev_i spend_ic\n' ...
             '  Then rerun Dynare (with bayesian_irf) and this script.'], ...
             strjoin(missing_own, ', '));
else
    fig = figure('Position', [100 100 900 550], 'Visible', 'off', ...
                 'Color', 'white', 'PaperPositionMode', 'auto');
    hold on;

    legend_handles = [];
    legend_names   = {};

    d1 = disc(1:Y1_QTRS);

    for s = 1:length(shocks)
        shock   = shocks{s};
        own_var = own_rev_map.(shock);

        y_fld   = ['y_' shock];
        own_fld = [own_var '_' shock];
        y_irf   = Median.(y_fld)(:);
        own_irf = Median.(own_fld)(:);

        % Year-1 (first 4 quarters) discounted cumulative own-revenue response
        den_y1 = sum(d1 .* own_irf(1:Y1_QTRS));

        if abs(den_y1) < 1e-12
            mult = NaN(T,1);
        else
            % Sign convention (matches compute_multipliers.m):
            %   Spending/ITC: no minus (expansionary shock, DY and Dspend same sign)
            %   Tax rates:    minus (positive shock raises tax, DY < 0)
            if strcmp(shock, 'eps_g') || strcmp(shock, 'eps_ic')
                mult = cumsum(disc .* y_irf) / den_y1;
            else
                mult = -cumsum(disc .* y_irf) / den_y1;
            end
        end

        h = plot(t, mult, '-', 'Color', line_colors(s,:), 'LineWidth', 2);
        legend_handles(end+1) = h; %#ok<SAGROW>
        legend_names{end+1}   = shock_labels{s}; %#ok<SAGROW>
    end

    plot([1 T], [0 0], '--', 'Color', zero_color, 'LineWidth', 0.5);
    hold off;

    xlabel('Horizon (quarters)', 'FontSize', 11);
    ylabel('Present-Value Multiplier', 'FontSize', 11);
    title('Fiscal Multipliers by Horizon (Year-1 Own-Revenue Denominator)', ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend(legend_handles, legend_names, 'Location', 'best', 'FontSize', 10);
    grid on;
    xlim([1 T]);
    set(gca, 'FontSize', 10, 'Box', 'on');

    filename = fullfile(fig_dir, 'multipliers_over_time_rev.png');
    print(fig, filename, '-dpng', '-r300');
    close(fig);
    fprintf('Saved %s\n', filename);
end

fprintf('\nAll figures written to %s\\ subfolder.\n', fig_dir);
