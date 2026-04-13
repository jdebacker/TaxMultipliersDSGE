% compute_multipliers_from_mcmc.m
% =========================================================================
% Computes present-value fiscal multipliers from posterior MCMC draws.
% Extended version: adds consumption tax, ITC, interest income tax,
% dividend tax, expensing rate, and transfer multipliers.
%
% Workflow:
%   1. Read posterior draws from the metropolis folder
%   2. Randomly select draws after burn-in
%   3. Write each draw into the model with set_all_parameters
%   4. Re-solve the model at that draw
%   5. Compute IRFs with simult_
%   6. Report median and [5,95] percentiles
%
% Usage:
%   dynare tax_dsge_est noclearall
%   compute_multipliers_from_mcmc
%
% Requires M_, oo_, options_, estim_params_ in the workspace.
% =========================================================================

if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ~exist('estim_params_','var')
    error('Run Dynare first: dynare tax_dsge_est noclearall');
end

%% ========================================================================
%  CONFIGURATION
%% ========================================================================
N_draws = 1000;
H       = [1 4 12 20];
T_irf   = max(H);
clip_val = 50;

if isfield(options_, 'mh_drop') && ~isempty(options_.mh_drop)
    mh_drop = options_.mh_drop;
else
    mh_drop = 0.333;
end

rng(2014);

%% ========================================================================
%  SAVE BASELINE STRUCTURES
%% ========================================================================
M_base       = M_;
oo_base      = oo_;
options_base = options_;

%% ========================================================================
%  LOAD POSTERIOR DRAWS
%% ========================================================================
fprintf('Loading MCMC draws from metropolis folder...\n');

mhpath = fullfile(M_.dname, 'metropolis');
if ~exist(mhpath, 'dir')
    error('Metropolis folder not found: %s', mhpath);
end

hist_file = fullfile(mhpath, [M_.fname '_mh_history_0.mat']);
if ~exist(hist_file, 'file')
    error('MH history file not found: %s', hist_file);
end
load(hist_file, 'record');

n_blocks = record.Nblck;
n_files  = record.LastFileNumber;

fprintf('  Found %d chain(s), %d file(s) per chain.\n', n_blocks, n_files);

post_draws = read_posterior_draws(mhpath, M_.fname, n_blocks, n_files, mh_drop);
if isempty(post_draws)
    error('No posterior draws loaded from %s', mhpath);
end

n_post = size(post_draws, 1);
n_cols = size(post_draws, 2);
fprintf('  Posterior draws kept after burn-in: %d\n', n_post);
fprintf('  Draw vector length: %d\n', n_cols);

if n_post < N_draws
    warning('Only %d posterior draws available; using all of them.', n_post);
    N_draws = n_post;
end

sel_idx   = randperm(n_post, N_draws);
sel_draws = post_draws(sel_idx, :);

fprintf('  Selected %d draws for multiplier computation.\n\n', N_draws);

%% ========================================================================
%  PRE-COMPUTE INDICES
%% ========================================================================
vnames = cellstr(M_.endo_names);
snames = cellstr(M_.exo_names);

% --- Shock indices ---
si_eps_g    = find(strcmp(snames, 'eps_g'));
si_eps_w    = find(strcmp(snames, 'eps_w'));
si_eps_k    = find(strcmp(snames, 'eps_k'));
si_eps_tr   = find(strcmp(snames, 'eps_tr'));
si_eps_c    = find(strcmp(snames, 'eps_c'));
si_eps_ic   = find(strcmp(snames, 'eps_ic'));
si_eps_ti   = find(strcmp(snames, 'eps_ti'));
si_eps_td   = find(strcmp(snames, 'eps_td'));
si_eps_etau = find(strcmp(snames, 'eps_etau'));

if isempty(si_eps_g) || isempty(si_eps_w) || isempty(si_eps_k)
    error('Could not find shock indices for eps_g, eps_w, eps_k.');
end
if isempty(si_eps_tr) || isempty(si_eps_c) || isempty(si_eps_ic) || ...
   isempty(si_eps_ti) || isempty(si_eps_td) || isempty(si_eps_etau)
    error('Could not find shock indices for extension shocks.');
end

% --- Endogenous variable indices ---
get_vi = @(nm) find(strcmp(vnames, nm));
vi_y      = get_vi('y');
vi_g      = get_vi('g');
vi_tr     = get_vi('tr');
vi_taxrev = get_vi('taxrev');
vi_R      = get_vi('R');

if any([isempty(vi_y), isempty(vi_g), isempty(vi_tr), isempty(vi_taxrev), isempty(vi_R)])
    error('Could not find required endogenous variables y, g, tr, taxrev, or R.');
end

if n_cols ~= number_of_estimated_objects(estim_params_)
    error(['MH draw length (%d) does not match Dynare estimated object count (%d). ' ...
           'Check the metropolis files or estimation setup.'], ...
           n_cols, number_of_estimated_objects(estim_params_));
end

%% ========================================================================
%  STORAGE
%% ========================================================================
% Zubairy original multipliers
multipliers_g  = NaN(N_draws, numel(H));   % govt spending
multipliers_w  = NaN(N_draws, numel(H));   % labor tax
multipliers_k  = NaN(N_draws, numel(H));   % capital tax

% Extension multipliers
multipliers_tr   = NaN(N_draws, numel(H)); % transfers
multipliers_c    = NaN(N_draws, numel(H)); % consumption tax
multipliers_ic   = NaN(N_draws, numel(H)); % investment tax credit
multipliers_ti   = NaN(N_draws, numel(H)); % interest income tax
multipliers_td   = NaN(N_draws, numel(H)); % dividend tax
multipliers_etau = NaN(N_draws, numel(H)); % expensing rate

n_success = 0;
n_fail    = 0;
t_start   = tic;

fprintf('Computing multipliers for %d posterior draws...\n', N_draws);

%% ========================================================================
%  MAIN LOOP
%% ========================================================================
for d = 1:N_draws
    try
        M_draw  = M_base;
        oo_draw = oo_base;

        xparam1 = sel_draws(d, :)';
        M_draw  = set_all_parameters(xparam1, estim_params_, M_draw);

        [dr_d, info, M_draw, oo_draw] = resol(0, M_draw, options_base, oo_draw);
        if info(1) ~= 0
            n_fail = n_fail + 1;
            continue;
        end

        ys_draw = oo_draw.dr.ys;
        R_ss    = ys_draw(vi_R);
        disc    = R_ss .^ (-(0:(T_irf-1))');

        % --- Zubairy original IRFs ---
        irf_g  = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_g, T_irf);
        irf_w  = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_w, T_irf);
        irf_k  = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_k, T_irf);

        % --- Extension IRFs ---
        irf_tr   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_tr, T_irf);
        irf_c    = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_c, T_irf);
        irf_ic   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_ic, T_irf);
        irf_ti   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_ti, T_irf);
        irf_td   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_td, T_irf);
        irf_etau = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_etau, T_irf);

        % --- Extract responses ---
        % Govt spending: PV(DY) / PV(DG)
        y_g   = irf_g(vi_y, :)';
        g_g   = irf_g(vi_g, :)';

        % Labor tax: -PV(DY) / PV(Dtaxrev)
        y_w   = irf_w(vi_y, :)';
        tr_w  = irf_w(vi_taxrev, :)';

        % Capital tax: -PV(DY) / PV(Dtaxrev)
        y_k   = irf_k(vi_y, :)';
        tr_k  = irf_k(vi_taxrev, :)';

        % Transfers: PV(DY) / PV(Dtr)
        y_tr    = irf_tr(vi_y, :)';
        tr_tr   = irf_tr(vi_tr, :)';

        % Consumption tax: -PV(DY) / PV(Dtaxrev)
        y_c     = irf_c(vi_y, :)';
        tr_c    = irf_c(vi_taxrev, :)';

        % ITC: -PV(DY) / PV(Dtaxrev)
        y_ic    = irf_ic(vi_y, :)';
        tr_ic   = irf_ic(vi_taxrev, :)';

        % Interest income tax: -PV(DY) / PV(Dtaxrev)
        y_ti    = irf_ti(vi_y, :)';
        tr_ti   = irf_ti(vi_taxrev, :)';

        % Dividend tax: -PV(DY) / PV(Dtaxrev)
        y_td    = irf_td(vi_y, :)';
        tr_td   = irf_td(vi_taxrev, :)';

        % Expensing rate: -PV(DY) / PV(Dtaxrev)
        y_etau  = irf_etau(vi_y, :)';
        tr_etau = irf_etau(vi_taxrev, :)';

        % --- Discounted cumulative sums ---
        pv_y_g  = cumsum(disc .* y_g);
        pv_g_g  = cumsum(disc .* g_g);
        pv_y_w  = cumsum(disc .* y_w);
        pv_tr_w = cumsum(disc .* tr_w);
        pv_y_k  = cumsum(disc .* y_k);
        pv_tr_k = cumsum(disc .* tr_k);

        pv_y_tr   = cumsum(disc .* y_tr);
        pv_tr_tr  = cumsum(disc .* tr_tr);
        pv_y_c    = cumsum(disc .* y_c);
        pv_tr_c   = cumsum(disc .* tr_c);
        pv_y_ic   = cumsum(disc .* y_ic);
        pv_tr_ic  = cumsum(disc .* tr_ic);
        pv_y_ti   = cumsum(disc .* y_ti);
        pv_tr_ti  = cumsum(disc .* tr_ti);
        pv_y_td   = cumsum(disc .* y_td);
        pv_tr_td  = cumsum(disc .* tr_td);
        pv_y_etau = cumsum(disc .* y_etau);
        pv_tr_etau = cumsum(disc .* tr_etau);

        % --- Multiplier ratios ---
        % Spending-type: PV(DY) / PV(DSpending)
        m_g  = safe_ratio(pv_y_g,   pv_g_g);
        m_tr = safe_ratio(pv_y_tr,  pv_tr_tr);

        % Tax-type: -PV(DY) / PV(Dtaxrev)
        m_w    = safe_ratio(-pv_y_w,    pv_tr_w);
        m_k    = safe_ratio(-pv_y_k,    pv_tr_k);
        m_c    = safe_ratio(-pv_y_c,    pv_tr_c);
        m_ic   = safe_ratio(-pv_y_ic,   pv_tr_ic);
        m_ti   = safe_ratio(-pv_y_ti,   pv_tr_ti);
        m_td   = safe_ratio(-pv_y_td,   pv_tr_td);
        m_etau = safe_ratio(-pv_y_etau, pv_tr_etau);

        % --- Store ---
        multipliers_g(d, :)    = m_g(H)';
        multipliers_w(d, :)    = m_w(H)';
        multipliers_k(d, :)    = m_k(H)';
        multipliers_tr(d, :)   = m_tr(H)';
        multipliers_c(d, :)    = m_c(H)';
        multipliers_ic(d, :)   = m_ic(H)';
        multipliers_ti(d, :)   = m_ti(H)';
        multipliers_td(d, :)   = m_td(H)';
        multipliers_etau(d, :) = m_etau(H)';

        n_success = n_success + 1;

    catch ME
        n_fail = n_fail + 1;
        if n_fail <= 3
            fprintf('  Draw %d failed: %s\n', d, ME.message);
            if ~isempty(ME.stack)
                fprintf('    In: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
            end
        end
        continue;
    end

    if mod(d, 100) == 0
        fprintf('  Draw %d / %d  (success: %d, fail: %d)\n', d, N_draws, n_success, n_fail);
    end
end

elapsed = toc(t_start);
fprintf('\nDone. %d successful, %d failed (%.1f sec)\n\n', n_success, n_fail, elapsed);

%% ========================================================================
%  RESTORE BASELINE WORKSPACE
%% ========================================================================
M_       = M_base;
oo_      = oo_base;
options_ = options_base;

%% ========================================================================
%  CLEAN AND SUMMARIZE
%% ========================================================================
valid_g    = clean_multiplier_matrix(multipliers_g, clip_val);
valid_w    = clean_multiplier_matrix(multipliers_w, clip_val);
valid_k    = clean_multiplier_matrix(multipliers_k, clip_val);
valid_tr   = clean_multiplier_matrix(multipliers_tr, clip_val);
valid_c    = clean_multiplier_matrix(multipliers_c, clip_val);
valid_ic   = clean_multiplier_matrix(multipliers_ic, clip_val);
valid_ti   = clean_multiplier_matrix(multipliers_ti, clip_val);
valid_td   = clean_multiplier_matrix(multipliers_td, clip_val);
valid_etau = clean_multiplier_matrix(multipliers_etau, clip_val);

if isempty(valid_g) || isempty(valid_w) || isempty(valid_k)
    error('No valid multiplier draws remain after cleaning (Zubairy originals).');
end

% --- Zubairy original ---
med_g = column_nanmedian(valid_g);  p05_g = column_prctile(valid_g, 5);  p95_g = column_prctile(valid_g, 95);
med_w = column_nanmedian(valid_w);  p05_w = column_prctile(valid_w, 5);  p95_w = column_prctile(valid_w, 95);
med_k = column_nanmedian(valid_k);  p05_k = column_prctile(valid_k, 5);  p95_k = column_prctile(valid_k, 95);

% --- Extension ---
med_tr   = column_nanmedian(valid_tr);    p05_tr   = column_prctile(valid_tr, 5);    p95_tr   = column_prctile(valid_tr, 95);
med_c    = column_nanmedian(valid_c);     p05_c    = column_prctile(valid_c, 5);     p95_c    = column_prctile(valid_c, 95);
med_ic   = column_nanmedian(valid_ic);    p05_ic   = column_prctile(valid_ic, 5);    p95_ic   = column_prctile(valid_ic, 95);
med_ti   = column_nanmedian(valid_ti);    p05_ti   = column_prctile(valid_ti, 5);    p95_ti   = column_prctile(valid_ti, 95);
med_td   = column_nanmedian(valid_td);    p05_td   = column_prctile(valid_td, 5);    p95_td   = column_prctile(valid_td, 95);
med_etau = column_nanmedian(valid_etau);  p05_etau = column_prctile(valid_etau, 5);  p95_etau = column_prctile(valid_etau, 95);

%% ========================================================================
%  ZUBAIRY (2014) TABLE 2 VALUES
%% ========================================================================
z_med_g = [1.07 1.06 0.90 0.72];
z_p05_g = [1.01 0.96 0.73 0.49];
z_p95_g = [1.13 1.17 1.07 0.93];

z_med_w = [0.13 0.32 0.68 0.85];
z_p05_w = [0.09 0.22 0.41 0.42];
z_p95_w = [0.18 0.45 1.09 1.58];

z_med_k = [0.34 0.43 0.52 0.46];
z_p05_k = [0.30 0.34 0.30 0.15];
z_p95_k = [0.37 0.51 0.73 0.81];

%% ========================================================================
%  PRINT TABLE
%% ========================================================================
fprintf('=======================================================================\n');
fprintf('TABLE 2. Present Value Multipliers (from MCMC Posterior Draws)\n');
fprintf('=======================================================================\n');
fprintf('Horizon (quarters):                   %5d  %5d  %5d  %5d\n\n', H);

fprintf('Panel A. Zubairy-Original Instruments (%d draws, %d successful)\n', N_draws, n_success);
print_multiplier_row('Govt spending',  med_g, p05_g, p95_g);
print_multiplier_row('Labor tax',      med_w, p05_w, p95_w);
print_multiplier_row('Capital tax',    med_k, p05_k, p95_k);
fprintf('\n');

fprintf('Panel B. Extension Instruments\n');
print_multiplier_row('Transfers',      med_tr,   p05_tr,   p95_tr);
print_multiplier_row('Consumption tax', med_c,   p05_c,    p95_c);
print_multiplier_row('ITC',            med_ic,   p05_ic,   p95_ic);
print_multiplier_row('Interest inc tax', med_ti, p05_ti,   p95_ti);
print_multiplier_row('Dividend tax',   med_td,   p05_td,   p95_td);
print_multiplier_row('Expensing rate', med_etau, p05_etau, p95_etau);
fprintf('\n');

fprintf('Panel C. Zubairy (2014) Table 2\n');
print_multiplier_row('Govt spending',  z_med_g, z_p05_g, z_p95_g);
print_multiplier_row('Labor tax',      z_med_w, z_p05_w, z_p95_w);
print_multiplier_row('Capital tax',    z_med_k, z_p05_k, z_p95_k);
fprintf('=======================================================================\n');

%% ========================================================================
%  SAVE CSV
%% ========================================================================
fmt_med = @(x) arrayfun(@(v) sprintf('%.2f', v), x, 'UniformOutput', false);
fmt_ci  = @(lo,hi) arrayfun(@(a,b) sprintf('[%.2f,%.2f]', a, b), lo, hi, 'UniformOutput', false);

rows = {};
% Panel A
rows(end+1,:) = [{'A','Govt spending','median'}, fmt_med(med_g)];
rows(end+1,:) = [{'A','Govt spending','[5,95]'}, fmt_ci(p05_g,p95_g)];
rows(end+1,:) = [{'A','Labor tax','median'}, fmt_med(med_w)];
rows(end+1,:) = [{'A','Labor tax','[5,95]'}, fmt_ci(p05_w,p95_w)];
rows(end+1,:) = [{'A','Capital tax','median'}, fmt_med(med_k)];
rows(end+1,:) = [{'A','Capital tax','[5,95]'}, fmt_ci(p05_k,p95_k)];
% Panel B
rows(end+1,:) = [{'B','Transfers','median'}, fmt_med(med_tr)];
rows(end+1,:) = [{'B','Transfers','[5,95]'}, fmt_ci(p05_tr,p95_tr)];
rows(end+1,:) = [{'B','Consumption tax','median'}, fmt_med(med_c)];
rows(end+1,:) = [{'B','Consumption tax','[5,95]'}, fmt_ci(p05_c,p95_c)];
rows(end+1,:) = [{'B','ITC','median'}, fmt_med(med_ic)];
rows(end+1,:) = [{'B','ITC','[5,95]'}, fmt_ci(p05_ic,p95_ic)];
rows(end+1,:) = [{'B','Interest inc tax','median'}, fmt_med(med_ti)];
rows(end+1,:) = [{'B','Interest inc tax','[5,95]'}, fmt_ci(p05_ti,p95_ti)];
rows(end+1,:) = [{'B','Dividend tax','median'}, fmt_med(med_td)];
rows(end+1,:) = [{'B','Dividend tax','[5,95]'}, fmt_ci(p05_td,p95_td)];
rows(end+1,:) = [{'B','Expensing rate','median'}, fmt_med(med_etau)];
rows(end+1,:) = [{'B','Expensing rate','[5,95]'}, fmt_ci(p05_etau,p95_etau)];
% Panel C
rows(end+1,:) = [{'C','Govt spending','median'}, fmt_med(z_med_g)];
rows(end+1,:) = [{'C','Govt spending','[5,95]'}, fmt_ci(z_p05_g,z_p95_g)];
rows(end+1,:) = [{'C','Labor tax','median'}, fmt_med(z_med_w)];
rows(end+1,:) = [{'C','Labor tax','[5,95]'}, fmt_ci(z_p05_w,z_p95_w)];
rows(end+1,:) = [{'C','Capital tax','median'}, fmt_med(z_med_k)];
rows(end+1,:) = [{'C','Capital tax','[5,95]'}, fmt_ci(z_p05_k,z_p95_k)];

T_out = cell2table(rows, 'VariableNames', ...
    {'Panel','Multiplier','Statistic','H1','H4','H12','H20'});
writetable(T_out, 'table2_multipliers_mcmc.csv', 'Delimiter', ',', 'QuoteStrings', true);
fprintf('Results saved to table2_multipliers_mcmc.csv\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================
function print_multiplier_row(name, med, p05, p95)
    fprintf('  %-22s (median):   %5.2f  %5.2f  %5.2f  %5.2f\n', name, med);
    fprintf('  %-22s [5,95]:    [%.2f,%.2f] [%.2f,%.2f] [%.2f,%.2f] [%.2f,%.2f]\n', ...
        name, p05(1),p95(1), p05(2),p95(2), p05(3),p95(3), p05(4),p95(4));
end

function post_draws = read_posterior_draws(mhpath, fname_root, n_blocks, n_files, mh_drop)
    post_draws = [];
    for b = 1:n_blocks
        chain_draws = [];
        for f = 1:n_files
            mh_file = fullfile(mhpath, sprintf('%s_mh%d_blck%d.mat', fname_root, f, b));
            if ~exist(mh_file, 'file')
                continue;
            end
            S = load(mh_file, 'x2');
            if ~isfield(S, 'x2') || isempty(S.x2)
                continue;
            end
            chain_draws = [chain_draws; S.x2]; %#ok<AGROW>
        end

        if isempty(chain_draws)
            warning('No draws loaded for block %d.', b);
            continue;
        end

        burn = ceil(mh_drop * size(chain_draws, 1));
        if burn >= size(chain_draws, 1)
            warning('All draws discarded by burn-in for block %d.', b);
            continue;
        end

        post_draws = [post_draws; chain_draws(burn+1:end, :)]; %#ok<AGROW>
    end
end

function n = number_of_estimated_objects(estim_params_)
    n = size(estim_params_.var_exo, 1) ...
      + size(estim_params_.var_endo, 1) ...
      + size(estim_params_.corrx, 1) ...
      + size(estim_params_.corrn, 1) ...
      + size(estim_params_.param_vals, 1);
end

function irf = compute_irf_simult(M_, options_, dr, ys, shock_idx, T)
    ex_ = zeros(T, M_.exo_nbr);
    ex_(1, shock_idx) = 1;
    y_path = simult_(M_, options_, ys, dr, ex_, 1);

    if size(y_path, 2) == T + 1
        irf = y_path(:, 2:end) - ys;
    elseif size(y_path, 2) == T
        irf = y_path - ys;
    else
        error('Unexpected output size from simult_.');
    end
end

function out = safe_ratio(num, den)
    out = NaN(size(num));
    idx = abs(den) > 1e-10;
    out(idx) = num(idx) ./ den(idx);
end

function X_clean = clean_multiplier_matrix(X, clip_val)
    X(abs(X) > clip_val) = NaN;
    X_clean = X(~all(isnan(X), 2), :);
end

function med = column_nanmedian(X)
    med = NaN(1, size(X, 2));
    for j = 1:size(X, 2)
        xj = X(:, j);
        xj = xj(~isnan(xj));
        if ~isempty(xj)
            med(j) = median(xj);
        end
    end
end

function q = column_prctile(X, p)
    q = NaN(1, size(X, 2));
    for j = 1:size(X, 2)
        xj = X(:, j);
        xj = xj(~isnan(xj));
        if ~isempty(xj)
            q(j) = prctile(xj, p);
        end
    end
end
