% compute_multipliers.m
% =========================================================================
% Computes fiscal multipliers from posterior MCMC draws using two
% normalization methods, side-by-side from the same draws.
%
% Method 1 (base):
%     For tax shocks:      -PV(DY) / PV(D total_taxrev)   through horizon H
%     For spending shocks:  PV(DY) / PV(D spending)
%   Denominator = present-value cumulative change in TOTAL tax revenue.
%   Includes behavioral feedback across all revenue components.
%
% Method 2 (revenue):
%     For tax shocks:      -PV(DY) / D(own_rev)_year1
%     For ITC:              PV(DY) / D(spend_ic)_year1
%     For spending shocks:  PV(DY) / D(spending)_year1
%   Denominator = year-1 (first 4 quarters, discounted) cumulative change
%   in the INSTRUMENT'S OWN revenue/spending series.
%   Approximates ex-ante first-year dollar cost of the specific policy.
%
% Instrument-to-own-revenue mapping (Method 2):
%   eps_g    -> g          (spending)
%   eps_tr   -> tr         (spending)
%   eps_w    -> rev_w      (labor tax revenue)
%   eps_k    -> rev_k      (capital tax revenue)
%   eps_c    -> rev_c      (consumption tax revenue)
%   eps_ic   -> spend_ic   (ITC fiscal cost)
%   eps_ti   -> rev_i      (interest income tax revenue)
%   eps_td   -> rev_d      (dividend tax revenue)
%   eps_etau -> rev_k      (expensing operates through capital tax base)
%
% Usage:
%   dynare tax_dsge_est_test1 noclearall
%   compute_multipliers
%
% Output files:
%   table2_multipliers_mcmc_base.csv    -- Method 1 (Panels A, B, C)
%   table2_multipliers_mcmc_rev.csv     -- Method 2 (Panels A, B)
% =========================================================================

if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ~exist('estim_params_','var')
    error('Run Dynare first: dynare tax_dsge_est_test1 noclearall');
end

%% ========================================================================
%  CONFIGURATION
%% ========================================================================
N_draws  = 1000;
H        = [1 4 12 20];
T_irf    = max(H);
clip_val = 50;
Y1_QTRS  = 4;   % first-year window for Method 2 denominator

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
vi_y        = get_vi('y');
vi_g        = get_vi('g');
vi_tr       = get_vi('tr');
vi_taxrev   = get_vi('taxrev');
vi_R        = get_vi('R');
vi_rev_c    = get_vi('rev_c');
vi_rev_w    = get_vi('rev_w');
vi_rev_k    = get_vi('rev_k');
vi_rev_d    = get_vi('rev_d');
vi_rev_i    = get_vi('rev_i');
vi_spend_ic = get_vi('spend_ic');

if any([isempty(vi_y), isempty(vi_g), isempty(vi_tr), isempty(vi_taxrev), isempty(vi_R)])
    error('Could not find required endogenous variables y, g, tr, taxrev, or R.');
end
if any([isempty(vi_rev_c), isempty(vi_rev_w), isempty(vi_rev_k), ...
        isempty(vi_rev_d), isempty(vi_rev_i), isempty(vi_spend_ic)])
    error('Could not find own-revenue variables (rev_c, rev_w, rev_k, rev_d, rev_i, spend_ic).');
end

if n_cols ~= number_of_estimated_objects(estim_params_)
    error('MH draw length (%d) does not match Dynare estimated object count (%d).', ...
           n_cols, number_of_estimated_objects(estim_params_));
end

%% ========================================================================
%  STORAGE
%% ========================================================================
% Method 1: total-taxrev denominator
m1_g    = NaN(N_draws, numel(H));
m1_w    = NaN(N_draws, numel(H));
m1_k    = NaN(N_draws, numel(H));
m1_tr   = NaN(N_draws, numel(H));
m1_c    = NaN(N_draws, numel(H));
m1_ic   = NaN(N_draws, numel(H));
m1_ti   = NaN(N_draws, numel(H));
m1_td   = NaN(N_draws, numel(H));
m1_etau = NaN(N_draws, numel(H));

% Method 2: year-1 own-revenue denominator
m2_g    = NaN(N_draws, numel(H));
m2_w    = NaN(N_draws, numel(H));
m2_k    = NaN(N_draws, numel(H));
m2_tr   = NaN(N_draws, numel(H));
m2_c    = NaN(N_draws, numel(H));
m2_ic   = NaN(N_draws, numel(H));
m2_ti   = NaN(N_draws, numel(H));
m2_td   = NaN(N_draws, numel(H));
m2_etau = NaN(N_draws, numel(H));

n_success = 0;
n_fail    = 0;
t_start   = tic;

fprintf('Computing multipliers for %d posterior draws (Methods 1 and 2)...\n', N_draws);

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

        % --- IRFs (one per shock) ---
        irf_g    = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_g, T_irf);
        irf_w    = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_w, T_irf);
        irf_k    = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_k, T_irf);
        irf_tr   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_tr, T_irf);
        irf_c    = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_c, T_irf);
        irf_ic   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_ic, T_irf);
        irf_ti   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_ti, T_irf);
        irf_td   = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_td, T_irf);
        irf_etau = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, si_eps_etau, T_irf);

        % --- PV(DY) cumulative through each horizon (shared numerator) ---
        pv_y_g    = cumsum(disc .* irf_g(vi_y, :)');
        pv_y_w    = cumsum(disc .* irf_w(vi_y, :)');
        pv_y_k    = cumsum(disc .* irf_k(vi_y, :)');
        pv_y_tr   = cumsum(disc .* irf_tr(vi_y, :)');
        pv_y_c    = cumsum(disc .* irf_c(vi_y, :)');
        pv_y_ic   = cumsum(disc .* irf_ic(vi_y, :)');
        pv_y_ti   = cumsum(disc .* irf_ti(vi_y, :)');
        pv_y_td   = cumsum(disc .* irf_td(vi_y, :)');
        pv_y_etau = cumsum(disc .* irf_etau(vi_y, :)');

        % ----------------------------------------------------------------
        % METHOD 1: denominator = PV(D total taxrev) [or spending] through H
        % ----------------------------------------------------------------
        pv_g_g      = cumsum(disc .* irf_g(vi_g, :)');
        pv_tr_tr    = cumsum(disc .* irf_tr(vi_tr, :)');
        pv_tax_w    = cumsum(disc .* irf_w(vi_taxrev, :)');
        pv_tax_k    = cumsum(disc .* irf_k(vi_taxrev, :)');
        pv_tax_c    = cumsum(disc .* irf_c(vi_taxrev, :)');
        pv_tax_ic   = cumsum(disc .* irf_ic(vi_taxrev, :)');
        pv_tax_ti   = cumsum(disc .* irf_ti(vi_taxrev, :)');
        pv_tax_td   = cumsum(disc .* irf_td(vi_taxrev, :)');
        pv_tax_etau = cumsum(disc .* irf_etau(vi_taxrev, :)');

        v1_g    = safe_ratio( pv_y_g,    pv_g_g);
        v1_tr   = safe_ratio( pv_y_tr,   pv_tr_tr);
        v1_w    = safe_ratio(-pv_y_w,    pv_tax_w);
        v1_k    = safe_ratio(-pv_y_k,    pv_tax_k);
        v1_c    = safe_ratio(-pv_y_c,    pv_tax_c);
        v1_ic   = safe_ratio(-pv_y_ic,   pv_tax_ic);
        v1_ti   = safe_ratio(-pv_y_ti,   pv_tax_ti);
        v1_td   = safe_ratio(-pv_y_td,   pv_tax_td);
        v1_etau = safe_ratio(-pv_y_etau, pv_tax_etau);

        % ----------------------------------------------------------------
        % METHOD 2: denominator = year-1 own-revenue (a single scalar)
        % ----------------------------------------------------------------
        d1 = disc(1:Y1_QTRS);

        own_g_y1    = sum(d1 .* irf_g(vi_g,         1:Y1_QTRS)');
        own_tr_y1   = sum(d1 .* irf_tr(vi_tr,       1:Y1_QTRS)');
        own_w_y1    = sum(d1 .* irf_w(vi_rev_w,     1:Y1_QTRS)');
        own_k_y1    = sum(d1 .* irf_k(vi_rev_k,     1:Y1_QTRS)');
        own_c_y1    = sum(d1 .* irf_c(vi_rev_c,     1:Y1_QTRS)');
        own_ic_y1   = sum(d1 .* irf_ic(vi_spend_ic, 1:Y1_QTRS)');
        own_ti_y1   = sum(d1 .* irf_ti(vi_rev_i,    1:Y1_QTRS)');
        own_td_y1   = sum(d1 .* irf_td(vi_rev_d,    1:Y1_QTRS)');
        own_etau_y1 = sum(d1 .* irf_etau(vi_rev_k,  1:Y1_QTRS)');

        % Sign convention:
        %   Spending instruments (g, tr, ITC): NO minus. DY and Dspend move
        %     in same direction for an expansionary shock, ratio positive.
        %   Tax-rate instruments (w, k, c, ti, td): minus. Positive shock
        %     raises tax rate and revenue but lowers DY; -DY/DRev > 0.
        %   Expensing: rev_k drops when eps_etau rises; DY > 0;
        %     -DY/(negative) = positive. Same formula as other tax rates.
        ones_T = ones(T_irf, 1);
        v2_g    = safe_ratio( pv_y_g,    own_g_y1    * ones_T);
        v2_tr   = safe_ratio( pv_y_tr,   own_tr_y1   * ones_T);
        v2_w    = safe_ratio(-pv_y_w,    own_w_y1    * ones_T);
        v2_k    = safe_ratio(-pv_y_k,    own_k_y1    * ones_T);
        v2_c    = safe_ratio(-pv_y_c,    own_c_y1    * ones_T);
        v2_ic   = safe_ratio( pv_y_ic,   own_ic_y1   * ones_T);  % NO MINUS
        v2_ti   = safe_ratio(-pv_y_ti,   own_ti_y1   * ones_T);
        v2_td   = safe_ratio(-pv_y_td,   own_td_y1   * ones_T);
        v2_etau = safe_ratio(-pv_y_etau, own_etau_y1 * ones_T);

        % --- Store at horizons ---
        m1_g(d,:)    = v1_g(H)';    m2_g(d,:)    = v2_g(H)';
        m1_w(d,:)    = v1_w(H)';    m2_w(d,:)    = v2_w(H)';
        m1_k(d,:)    = v1_k(H)';    m2_k(d,:)    = v2_k(H)';
        m1_tr(d,:)   = v1_tr(H)';   m2_tr(d,:)   = v2_tr(H)';
        m1_c(d,:)    = v1_c(H)';    m2_c(d,:)    = v2_c(H)';
        m1_ic(d,:)   = v1_ic(H)';   m2_ic(d,:)   = v2_ic(H)';
        m1_ti(d,:)   = v1_ti(H)';   m2_ti(d,:)   = v2_ti(H)';
        m1_td(d,:)   = v1_td(H)';   m2_td(d,:)   = v2_td(H)';
        m1_etau(d,:) = v1_etau(H)'; m2_etau(d,:) = v2_etau(H)';

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
% Method 1
valid_m1_g    = clean_multiplier_matrix(m1_g,    clip_val);
valid_m1_w    = clean_multiplier_matrix(m1_w,    clip_val);
valid_m1_k    = clean_multiplier_matrix(m1_k,    clip_val);
valid_m1_tr   = clean_multiplier_matrix(m1_tr,   clip_val);
valid_m1_c    = clean_multiplier_matrix(m1_c,    clip_val);
valid_m1_ic   = clean_multiplier_matrix(m1_ic,   clip_val);
valid_m1_ti   = clean_multiplier_matrix(m1_ti,   clip_val);
valid_m1_td   = clean_multiplier_matrix(m1_td,   clip_val);
valid_m1_etau = clean_multiplier_matrix(m1_etau, clip_val);

% Method 2
valid_m2_g    = clean_multiplier_matrix(m2_g,    clip_val);
valid_m2_w    = clean_multiplier_matrix(m2_w,    clip_val);
valid_m2_k    = clean_multiplier_matrix(m2_k,    clip_val);
valid_m2_tr   = clean_multiplier_matrix(m2_tr,   clip_val);
valid_m2_c    = clean_multiplier_matrix(m2_c,    clip_val);
valid_m2_ic   = clean_multiplier_matrix(m2_ic,   clip_val);
valid_m2_ti   = clean_multiplier_matrix(m2_ti,   clip_val);
valid_m2_td   = clean_multiplier_matrix(m2_td,   clip_val);
valid_m2_etau = clean_multiplier_matrix(m2_etau, clip_val);

if isempty(valid_m1_g) || isempty(valid_m1_w) || isempty(valid_m1_k)
    error('No valid multiplier draws remain after cleaning.');
end

% Method 1 summaries
med1_g = column_nanmedian(valid_m1_g);  p05_1_g = column_prctile(valid_m1_g, 5);  p95_1_g = column_prctile(valid_m1_g, 95);
med1_w = column_nanmedian(valid_m1_w);  p05_1_w = column_prctile(valid_m1_w, 5);  p95_1_w = column_prctile(valid_m1_w, 95);
med1_k = column_nanmedian(valid_m1_k);  p05_1_k = column_prctile(valid_m1_k, 5);  p95_1_k = column_prctile(valid_m1_k, 95);
med1_tr   = column_nanmedian(valid_m1_tr);    p05_1_tr   = column_prctile(valid_m1_tr, 5);    p95_1_tr   = column_prctile(valid_m1_tr, 95);
med1_c    = column_nanmedian(valid_m1_c);     p05_1_c    = column_prctile(valid_m1_c, 5);     p95_1_c    = column_prctile(valid_m1_c, 95);
med1_ic   = column_nanmedian(valid_m1_ic);    p05_1_ic   = column_prctile(valid_m1_ic, 5);    p95_1_ic   = column_prctile(valid_m1_ic, 95);
med1_ti   = column_nanmedian(valid_m1_ti);    p05_1_ti   = column_prctile(valid_m1_ti, 5);    p95_1_ti   = column_prctile(valid_m1_ti, 95);
med1_td   = column_nanmedian(valid_m1_td);    p05_1_td   = column_prctile(valid_m1_td, 5);    p95_1_td   = column_prctile(valid_m1_td, 95);
med1_etau = column_nanmedian(valid_m1_etau);  p05_1_etau = column_prctile(valid_m1_etau, 5);  p95_1_etau = column_prctile(valid_m1_etau, 95);

% Method 2 summaries
med2_g = column_nanmedian(valid_m2_g);  p05_2_g = column_prctile(valid_m2_g, 5);  p95_2_g = column_prctile(valid_m2_g, 95);
med2_w = column_nanmedian(valid_m2_w);  p05_2_w = column_prctile(valid_m2_w, 5);  p95_2_w = column_prctile(valid_m2_w, 95);
med2_k = column_nanmedian(valid_m2_k);  p05_2_k = column_prctile(valid_m2_k, 5);  p95_2_k = column_prctile(valid_m2_k, 95);
med2_tr   = column_nanmedian(valid_m2_tr);    p05_2_tr   = column_prctile(valid_m2_tr, 5);    p95_2_tr   = column_prctile(valid_m2_tr, 95);
med2_c    = column_nanmedian(valid_m2_c);     p05_2_c    = column_prctile(valid_m2_c, 5);     p95_2_c    = column_prctile(valid_m2_c, 95);
med2_ic   = column_nanmedian(valid_m2_ic);    p05_2_ic   = column_prctile(valid_m2_ic, 5);    p95_2_ic   = column_prctile(valid_m2_ic, 95);
med2_ti   = column_nanmedian(valid_m2_ti);    p05_2_ti   = column_prctile(valid_m2_ti, 5);    p95_2_ti   = column_prctile(valid_m2_ti, 95);
med2_td   = column_nanmedian(valid_m2_td);    p05_2_td   = column_prctile(valid_m2_td, 5);    p95_2_td   = column_prctile(valid_m2_td, 95);
med2_etau = column_nanmedian(valid_m2_etau);  p05_2_etau = column_prctile(valid_m2_etau, 5);  p95_2_etau = column_prctile(valid_m2_etau, 95);

%% ========================================================================
%  ZUBAIRY (2014) TABLE 2 REFERENCE VALUES  (Panel C, Method 1 CSV only)
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
%  PRINT TABLES
%% ========================================================================
% --- Method 1 ---
fprintf('=======================================================================\n');
fprintf('TABLE 2. Method 1: -PV(DY) / PV(D total taxrev)\n');
fprintf('=======================================================================\n');
fprintf('Horizon (quarters):                   %5d  %5d  %5d  %5d\n\n', H);
fprintf('Panel A. Zubairy-Original Instruments (%d draws, %d successful)\n', N_draws, n_success);
print_multiplier_row('Govt spending',  med1_g, p05_1_g, p95_1_g);
print_multiplier_row('Labor tax',      med1_w, p05_1_w, p95_1_w);
print_multiplier_row('Capital tax',    med1_k, p05_1_k, p95_1_k);
fprintf('\nPanel B. Extension Instruments\n');
print_multiplier_row('Transfers',       med1_tr,   p05_1_tr,   p95_1_tr);
print_multiplier_row('Consumption tax', med1_c,    p05_1_c,    p95_1_c);
print_multiplier_row('ITC',             med1_ic,   p05_1_ic,   p95_1_ic);
print_multiplier_row('Interest inc tax',med1_ti,   p05_1_ti,   p95_1_ti);
print_multiplier_row('Dividend tax',    med1_td,   p05_1_td,   p95_1_td);
print_multiplier_row('Expensing rate',  med1_etau, p05_1_etau, p95_1_etau);
fprintf('\nPanel C. Zubairy (2014) Table 2\n');
print_multiplier_row('Govt spending',  z_med_g, z_p05_g, z_p95_g);
print_multiplier_row('Labor tax',      z_med_w, z_p05_w, z_p95_w);
print_multiplier_row('Capital tax',    z_med_k, z_p05_k, z_p95_k);
fprintf('=======================================================================\n\n');

% --- Method 2 ---
fprintf('=======================================================================\n');
fprintf('TABLE 2. Method 2: -PV(DY) / D(own_rev)_year1\n');
fprintf('=======================================================================\n');
fprintf('Horizon (quarters):                   %5d  %5d  %5d  %5d\n\n', H);
fprintf('Panel A. Zubairy-Original Instruments\n');
print_multiplier_row('Govt spending',  med2_g, p05_2_g, p95_2_g);
print_multiplier_row('Labor tax',      med2_w, p05_2_w, p95_2_w);
print_multiplier_row('Capital tax',    med2_k, p05_2_k, p95_2_k);
fprintf('\nPanel B. Extension Instruments\n');
print_multiplier_row('Transfers',       med2_tr,   p05_2_tr,   p95_2_tr);
print_multiplier_row('Consumption tax', med2_c,    p05_2_c,    p95_2_c);
print_multiplier_row('ITC',             med2_ic,   p05_2_ic,   p95_2_ic);
print_multiplier_row('Interest inc tax',med2_ti,   p05_2_ti,   p95_2_ti);
print_multiplier_row('Dividend tax',    med2_td,   p05_2_td,   p95_2_td);
print_multiplier_row('Expensing rate',  med2_etau, p05_2_etau, p95_2_etau);
fprintf('=======================================================================\n\n');

%% ========================================================================
%  SAVE CSVs
%% ========================================================================
fmt_med = @(x) arrayfun(@(v) sprintf('%.2f', v), x, 'UniformOutput', false);
fmt_ci  = @(lo,hi) arrayfun(@(a,b) sprintf('[%.2f,%.2f]', a, b), lo, hi, 'UniformOutput', false);

% --- Method 1 CSV: Panels A, B, C ---
rows1 = {};
rows1(end+1,:) = [{'A','Govt spending','median'},    fmt_med(med1_g)];
rows1(end+1,:) = [{'A','Govt spending','[5,95]'},    fmt_ci(p05_1_g, p95_1_g)];
rows1(end+1,:) = [{'A','Labor tax','median'},        fmt_med(med1_w)];
rows1(end+1,:) = [{'A','Labor tax','[5,95]'},        fmt_ci(p05_1_w, p95_1_w)];
rows1(end+1,:) = [{'A','Capital tax','median'},      fmt_med(med1_k)];
rows1(end+1,:) = [{'A','Capital tax','[5,95]'},      fmt_ci(p05_1_k, p95_1_k)];
rows1(end+1,:) = [{'B','Transfers','median'},        fmt_med(med1_tr)];
rows1(end+1,:) = [{'B','Transfers','[5,95]'},        fmt_ci(p05_1_tr, p95_1_tr)];
rows1(end+1,:) = [{'B','Consumption tax','median'},  fmt_med(med1_c)];
rows1(end+1,:) = [{'B','Consumption tax','[5,95]'},  fmt_ci(p05_1_c, p95_1_c)];
rows1(end+1,:) = [{'B','ITC','median'},              fmt_med(med1_ic)];
rows1(end+1,:) = [{'B','ITC','[5,95]'},              fmt_ci(p05_1_ic, p95_1_ic)];
rows1(end+1,:) = [{'B','Interest inc tax','median'}, fmt_med(med1_ti)];
rows1(end+1,:) = [{'B','Interest inc tax','[5,95]'}, fmt_ci(p05_1_ti, p95_1_ti)];
rows1(end+1,:) = [{'B','Dividend tax','median'},     fmt_med(med1_td)];
rows1(end+1,:) = [{'B','Dividend tax','[5,95]'},     fmt_ci(p05_1_td, p95_1_td)];
rows1(end+1,:) = [{'B','Expensing rate','median'},   fmt_med(med1_etau)];
rows1(end+1,:) = [{'B','Expensing rate','[5,95]'},   fmt_ci(p05_1_etau, p95_1_etau)];
rows1(end+1,:) = [{'C','Govt spending','median'},    fmt_med(z_med_g)];
rows1(end+1,:) = [{'C','Govt spending','[5,95]'},    fmt_ci(z_p05_g, z_p95_g)];
rows1(end+1,:) = [{'C','Labor tax','median'},        fmt_med(z_med_w)];
rows1(end+1,:) = [{'C','Labor tax','[5,95]'},        fmt_ci(z_p05_w, z_p95_w)];
rows1(end+1,:) = [{'C','Capital tax','median'},      fmt_med(z_med_k)];
rows1(end+1,:) = [{'C','Capital tax','[5,95]'},      fmt_ci(z_p05_k, z_p95_k)];

T1 = cell2table(rows1, 'VariableNames', {'Panel','Multiplier','Statistic','H1','H4','H12','H20'});
writetable(T1, 'table2_multipliers_mcmc_base.csv', 'Delimiter', ',', 'QuoteStrings', true);

% --- Method 2 CSV: Panels A, B ---
rows2 = {};
rows2(end+1,:) = [{'A','Govt spending','median'},    fmt_med(med2_g)];
rows2(end+1,:) = [{'A','Govt spending','[5,95]'},    fmt_ci(p05_2_g, p95_2_g)];
rows2(end+1,:) = [{'A','Labor tax','median'},        fmt_med(med2_w)];
rows2(end+1,:) = [{'A','Labor tax','[5,95]'},        fmt_ci(p05_2_w, p95_2_w)];
rows2(end+1,:) = [{'A','Capital tax','median'},      fmt_med(med2_k)];
rows2(end+1,:) = [{'A','Capital tax','[5,95]'},      fmt_ci(p05_2_k, p95_2_k)];
rows2(end+1,:) = [{'B','Transfers','median'},        fmt_med(med2_tr)];
rows2(end+1,:) = [{'B','Transfers','[5,95]'},        fmt_ci(p05_2_tr, p95_2_tr)];
rows2(end+1,:) = [{'B','Consumption tax','median'},  fmt_med(med2_c)];
rows2(end+1,:) = [{'B','Consumption tax','[5,95]'},  fmt_ci(p05_2_c, p95_2_c)];
rows2(end+1,:) = [{'B','ITC','median'},              fmt_med(med2_ic)];
rows2(end+1,:) = [{'B','ITC','[5,95]'},              fmt_ci(p05_2_ic, p95_2_ic)];
rows2(end+1,:) = [{'B','Interest inc tax','median'}, fmt_med(med2_ti)];
rows2(end+1,:) = [{'B','Interest inc tax','[5,95]'}, fmt_ci(p05_2_ti, p95_2_ti)];
rows2(end+1,:) = [{'B','Dividend tax','median'},     fmt_med(med2_td)];
rows2(end+1,:) = [{'B','Dividend tax','[5,95]'},     fmt_ci(p05_2_td, p95_2_td)];
rows2(end+1,:) = [{'B','Expensing rate','median'},   fmt_med(med2_etau)];
rows2(end+1,:) = [{'B','Expensing rate','[5,95]'},   fmt_ci(p05_2_etau, p95_2_etau)];

T2 = cell2table(rows2, 'VariableNames', {'Panel','Multiplier','Statistic','H1','H4','H12','H20'});
writetable(T2, 'table2_multipliers_mcmc_rev.csv', 'Delimiter', ',', 'QuoteStrings', true);

fprintf('Saved:\n');
fprintf('  table2_multipliers_mcmc_base.csv    (Method 1)\n');
fprintf('  table2_multipliers_mcmc_rev.csv     (Method 2)\n');

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
