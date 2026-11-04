% compute_multipliers.m
% =========================================================================
% Fiscal multipliers from posterior MCMC draws, under three normalizations
% computed side by side from the same draws.
%
% -------------------------------------------------------------------------
% MAIN: PRESENT-VALUE OWN-REVENUE MULTIPLIER
% -------------------------------------------------------------------------
%     tax instruments:      -PV(dY, 0..H) / PV(d own_revenue, 0..H)
%     spending, ITC, expensing:      PV(dY, 0..H) / PV(d own_outlay,  0..H)
%
%   Numerator and denominator are cumulated to the SAME horizon H, and the
%   denominator is the instrument's own revenue or outlay. This is the
%   standard present-value convention in the fiscal multiplier literature
%   and is the object reported in the main table.
%   The main table keeps the finite ratio for every draw. A separate table
%   reports, for each tax instrument and horizon, the share of draws that
%   are self-financing: positive expansionary GDP with nonpositive total
%   net revenue cost. table_multipliers_pvown_inf.csv applies the +Inf
%   reporting override to those draws instead. Spending and transfers are
%   not classified. The mechanical and total-revenue tables keep the +Inf
%   override.
%
% -------------------------------------------------------------------------
% SECONDARY: MECHANICAL OWN-COST NORMALIZATION
% -------------------------------------------------------------------------
%     signed PV(dY, 0..H) / PV(mechanical own cost, 0..H)
%
%   The denominator holds the tax base at its steady state and lets only
%   the statutory rate path move, so it is the no-behavior fiscal cost of
%   the instrument. This is the convention used in the companion
%   heterogeneous-agent paper, so this table is the one directly
%   comparable across the two papers.
%
%   The rate path per unit innovation is rho^(t-1) under each instrument's
%   own fiscal rule. Log rules (tau_k, tau_w, tau_c, tau_i, tau_d, g, tr)
%   give d(revenue)_t = revenue_ss * rho^(t-1). Level rules (tau_ic,
%   e_tau) give ITC outlays I_ss * rho^(t-1) and expensing deductions
%   tauk_ss * I_ss * rho^(t-1). Expensing cost also subtracts the later
%   revenue recovery from the smaller tax basis, at rate delta_tau.
%   The mechanical rate path omits endogenous fiscal feedback; the main
%   expensing denominator uses the simulated e_tau path instead.
%
% -------------------------------------------------------------------------
% ROBUSTNESS: TOTAL-REVENUE NORMALIZATION  (previously the main method)
% -------------------------------------------------------------------------
%     tax instruments:      -PV(dY, 0..H) / PV(d total taxrev, 0..H)
%     spending:              PV(dY, 0..H) / PV(d own_outlay,   0..H)
%
%   The denominator is the change in TOTAL tax revenue, so it nets the
%   direct revenue loss against revenue induced elsewhere in the tax
%   system. For an instrument close to self-financing the denominator
%   approaches zero and the ratio becomes unstable. Retained for
%   comparison with the previous draft and with Zubairy (2014).
%
% -------------------------------------------------------------------------
% Instrument to own-revenue mapping
%   eps_g    -> g          (outlay)
%   eps_tr   -> tr         (outlay)
%   eps_w    -> rev_w
%   eps_k    -> rev_k
%   eps_c    -> rev_c
%   eps_ic   -> spend_ic   (outlay)
%   eps_ti   -> rev_i
%   eps_td   -> rev_d
%   eps_etau -> spend_etau (incremental expensing deductions less revenue
%                            recovered from the resulting smaller tax basis)
%
% Usage:
%   dynare tax_dsge_est noclearall
%   compute_multipliers
%
% Output files:
%   table_multipliers_pvown.csv    -- MAIN (finite ratios)
%   table_self_financing.csv       -- share of self-financing draws
%   table_multipliers_pvown_inf.csv -- MAIN with the +Inf override
%   table_multipliers_mech.csv     -- secondary, matches the HANK paper
%   table_multipliers_totalrev.csv -- robustness (includes Zubairy Panel C)
% =========================================================================

% Shared reporting helpers live in the parent Codes directory.
addpath(fileparts(fileparts(mfilename('fullpath'))));

if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ~exist('estim_params_','var')
    error('Run Dynare first: dynare tax_dsge_est noclearall');
end

%% ========================================================================
%  CONFIGURATION
%% ========================================================================
N_draws  = 1000;
H        = [1 4 12 20];
T_irf    = max(H);
clip_val = 50;

% -------------------------------------------------------------------------
% ISOLATE ONE INSTRUMENT AT A TIME
%
% phi is the coefficient on eps_k in the labor tax rule:
%     log(tau_w/tauw_ss) = ... + phi*eps_k + eps_w
% At its posterior value a capital tax innovation also moves the labor tax
% rate on impact, so the capital tax experiment is not a single-instrument
% experiment. Setting phi to zero here removes that spillover.
%
% This is a SIMULATION-ONLY override. phi stays in the model and in
% estimated_params, so the posterior is unchanged and the existing MCMC
% chains remain valid. No re-estimation is required.
%
% phi enters no other equation and does not enter the steady state, so this
% affects the eps_k experiment ONLY. Every other instrument should be
% numerically unchanged relative to a run with PHI_ZERO = false, which is a
% useful check that the override does what it should and nothing more.
% -------------------------------------------------------------------------
PHI_ZERO = true;

pnames_all = cellstr(M_.param_names);
ip_phi     = find(strcmp(pnames_all, 'phi'));
if PHI_ZERO && isempty(ip_phi)
    error('PHI_ZERO is set but parameter phi was not found in M_.param_names.');
end

if isfield(options_, 'mh_drop') && ~isempty(options_.mh_drop)
    mh_drop = options_.mh_drop;
else
    mh_drop = 0.333;
end

rng(2014);

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

post_draws = read_posterior_draws(mhpath, M_.fname, record.Nblck, ...
                                  record.LastFileNumber, mh_drop);
if isempty(post_draws)
    error('No posterior draws loaded from %s', mhpath);
end

n_post = size(post_draws, 1);
n_cols = size(post_draws, 2);
fprintf('  Posterior draws kept after burn-in: %d\n', n_post);

if n_cols ~= number_of_estimated_objects(estim_params_)
    error('MH draw length (%d) does not match Dynare estimated object count (%d).', ...
           n_cols, number_of_estimated_objects(estim_params_));
end

if n_post < N_draws
    warning('Only %d posterior draws available; using all of them.', n_post);
    N_draws = n_post;
end

sel_draws = post_draws(randperm(n_post, N_draws), :);
fprintf('  Selected %d draws.\n\n', N_draws);

%% ========================================================================
%  INDICES AND INSTRUMENT TABLE
%% ========================================================================
vnames = cellstr(M_.endo_names);
snames = cellstr(M_.exo_names);

gv = @(nm) idx_or_die(vnames, nm, 'endogenous variable');
gs = @(nm) idx_or_die(snames, nm, 'shock');

vi_y      = gv('y');
vi_R      = gv('R');
vi_taxrev = gv('taxrev');

% -------------------------------------------------------------------------
%   panel   : 'A' = Zubairy-original instruments, 'B' = extension
%   ownvar  : the instrument's own revenue or outlay series
%   numsign : sign applied to PV(dY) so an expansionary action gives a
%             positive multiplier.
%               +1 for outlays (g, tr, ITC, expensing): same sign as dY
%               -1 for tax rates: a rate increase raises revenue, lowers Y
%   totsign : sign applied to PV(dY) in the TOTAL-REVENUE table. Usually
%             equal to numsign, but ITC and expensing differ: their own outlay
%             RISES with the instrument while total revenue
%             FALLS, so the two denominators need opposite numerator signs.
%   usetax  : true  -> total-revenue robustness uses PV(d taxrev)
%             false -> total-revenue robustness uses the own outlay, since
%                      spending instruments have no meaningful total-
%                      revenue denominator
%   rho     : the instrument's own persistence parameter, used to build the
%             mechanical rate path
%   rule    : 'log'   -> d(own)_t = own_ss * rho^(t-1)
%             'ic'    -> d(spend_ic)_t = I_ss * rho^(t-1)
%             'etau'  -> current expensing cost less depreciation recovery
% -------------------------------------------------------------------------
INS = struct('label',{},'panel',{},'shock',{},'ownvar',{},'numsign',{}, ...
             'totsign',{},'usetax',{},'rho',{},'rule',{});
%                label              pnl shock      ownvar     num tot usetax rho        rule
INS(end+1) = mk('Govt spending',    'A','eps_g',   'g',        +1, +1,false,'rho_g',   'log' );
INS(end+1) = mk('Labor tax',        'A','eps_w',   'rev_w',    -1, -1,true ,'rho_w',   'log' );
INS(end+1) = mk('Capital tax',      'A','eps_k',   'rev_k',    -1, -1,true ,'rho_k',   'log' );
INS(end+1) = mk('Transfers',        'B','eps_tr',  'tr',       +1, +1,false,'rho_tr',  'log' );
INS(end+1) = mk('Consumption tax',  'B','eps_c',   'rev_c',    -1, -1,true ,'rho_c',   'log' );
INS(end+1) = mk('ITC',              'B','eps_ic',  'spend_ic', +1, -1,true ,'rho_ic',  'ic'  );
% Interest income tax exogenous and constant from 2026-10-22 (no eps_ti):
% INS(end+1) = mk('Interest inc tax', 'B','eps_ti',  'rev_i',    -1, -1,true ,'rho_ti',  'log' );
INS(end+1) = mk('Dividend tax',     'B','eps_td',  'rev_d',    -1, -1,true ,'rho_td',  'log' );
INS(end+1) = mk('Expensing rate',   'B','eps_etau','spend_etau',+1, -1,true ,'rho_etau','etau');

nI = numel(INS);
for j = 1:nI
    INS(j).si = gs(INS(j).shock);
    INS(j).vi = gv(INS(j).ownvar);
    INS(j).ir = idx_or_die(pnames_all, INS(j).rho, 'parameter');
end

vi_i       = gv('i');
ip_tauk_ss = idx_or_die(pnames_all, 'tauk_ss', 'parameter');
ip_delta_tau = idx_or_die(pnames_all, 'delta_tau', 'parameter');

%% ========================================================================
%  STORAGE
%% ========================================================================
m_pvown = NaN(N_draws, nI, numel(H));   % MAIN (finite ratios)
m_pvinf = NaN(N_draws, nI, numel(H));   % MAIN with the +Inf override
m_sf    = NaN(N_draws, nI, numel(H));   % 1 if self-financing, 0 if not
m_mech  = NaN(N_draws, nI, numel(H));   % secondary, matches HANK paper
m_total = NaN(N_draws, nI, numel(H));   % robustness

n_success = 0;
n_fail    = 0;
t_start   = tic;

fprintf('Computing multipliers over %d posterior draws...\n', N_draws);
if PHI_ZERO
    fprintf('  phi set to zero in the counterfactual (one instrument at a time).\n');
end

%% ========================================================================
%  MAIN LOOP
%% ========================================================================
for d = 1:N_draws
    try
        M_draw  = M_base;
        oo_draw = oo_base;

        M_draw = set_all_parameters(sel_draws(d,:)', estim_params_, M_draw);

        if PHI_ZERO
            M_draw.params(ip_phi) = 0;
        end

        [dr_d, info, M_draw, oo_draw] = resol(0, M_draw, options_base, oo_draw);
        if info(1) ~= 0
            n_fail = n_fail + 1;
            continue;
        end

        ys_draw = oo_draw.dr.ys;
        disc    = ys_draw(vi_R) .^ (-(0:(T_irf-1))');

        for j = 1:nI
            ins = INS(j);

            irf = compute_irf_simult(M_draw, options_base, dr_d, ys_draw, ins.si, T_irf);

            % Numerator, signed, cumulated to each horizon
            pv_y = ins.numsign * cumsum(disc .* irf(vi_y, :)');

            % --- MAIN: own revenue cumulated to the SAME horizon ---------
            pv_own = cumsum(disc .* irf(ins.vi, :)');
            % Classify using TOTAL revenue, regardless of normalization.
            % numsign also gives the expansionary shock direction: tax cuts
            % reverse rate innovations; ITC and expensing increases do not.
            net_cost = NaN(size(pv_y));
            if ins.usetax
                net_cost = -ins.numsign * cumsum(disc .* irf(vi_taxrev, :)');
            end
            v_pv     = safe_ratio(pv_y, pv_own, NaN(size(pv_y)));   % finite
            v_pv_inf = safe_ratio(pv_y, pv_own, net_cost);           % +Inf override
            if ins.usetax
                sf = double(isfinite(pv_y) & isfinite(net_cost) & pv_y > 1e-10 & net_cost <= 1e-10);
            else
                sf = NaN(size(pv_y));
            end

            % --- ROBUSTNESS: total revenue (own outlay for spending) -----
            % ITC and expensing need the opposite numerator sign: their cost
            % rises with the instrument while total revenue falls.
            if ins.usetax
                pv_den_tot = cumsum(disc .* irf(vi_taxrev, :)');
            else
                pv_den_tot = pv_own;
            end
            % Orient both inputs as expansionary GDP gain and fiscal cost.
            v_tot = safe_ratio(pv_y, ...
                (ins.totsign / ins.numsign) * pv_den_tot, net_cost);

            % --- SECONDARY: mechanical own cost --------------------------
            % Base frozen at steady state; only the statutory rate moves.
            rho_j = M_draw.params(ins.ir);
            decay = rho_j .^ (0:(T_irf-1))';
            switch ins.rule
                case 'log'
                    mech_path = ys_draw(ins.vi) * decay;
                case 'ic'
                    mech_path = ys_draw(vi_i) * decay;
                case 'etau'
                    % Tax-basis difference relative to unchanged expensing.
                    % Lagged basis is negative after an expensing increase,
                    % so its depreciation term reduces the net fiscal cost.
                    delta_tax = M_draw.params(ip_delta_tau);
                    basis_gap = 0;
                    mech_path = zeros(T_irf, 1);
                    for t = 1:T_irf
                        deduction = ys_draw(vi_i) * decay(t);
                        mech_path(t) = M_draw.params(ip_tauk_ss) ...
                            * (deduction + delta_tax*basis_gap);
                        basis_gap = (1 - delta_tax)*basis_gap - deduction;
                    end
            end
            pv_mech = cumsum(disc .* mech_path);
            v_mech  = safe_ratio(pv_y, pv_mech, net_cost);

            m_pvown(d,j,:) = v_pv(H)';
            m_pvinf(d,j,:) = v_pv_inf(H)';
            m_sf(d,j,:)    = sf(H)';
            m_mech(d,j,:)  = v_mech(H)';
            m_total(d,j,:) = v_tot(H)';
        end

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

fprintf('\nDone. %d successful, %d failed (%.1f sec)\n\n', n_success, n_fail, toc(t_start));

M_ = M_base; oo_ = oo_base; options_ = options_base;

if n_success == 0
    error('No draws solved successfully.');
end

%% ========================================================================
%  SUMMARIZE
%% ========================================================================
[med_pv, lo_pv, hi_pv] = summarize_block(m_pvown, clip_val);
[med_pi, lo_pi, hi_pi] = summarize_block(m_pvinf, clip_val);
sf_share = reshape(mean(m_sf, 1, 'omitnan'), nI, numel(H));   % NaN for spending, transfers
[med_mc, lo_mc, hi_mc] = summarize_block(m_mech,  clip_val);
[med_tt, lo_tt, hi_tt] = summarize_block(m_total, clip_val);

%% ========================================================================
%  ZUBAIRY (2014) TABLE 2 REFERENCE VALUES
%  Hardcoded reference numbers from the published paper, not model output.
%  Reported as Panel C of the total-revenue table for comparability.
%% ========================================================================
z_med_g = [1.07 1.06 0.90 0.72];  z_p05_g = [1.01 0.96 0.73 0.49];  z_p95_g = [1.13 1.17 1.07 0.93];
z_med_w = [0.13 0.32 0.68 0.85];  z_p05_w = [0.09 0.22 0.41 0.42];  z_p95_w = [0.18 0.45 1.09 1.58];
z_med_k = [0.34 0.43 0.52 0.46];  z_p05_k = [0.30 0.34 0.30 0.15];  z_p95_k = [0.37 0.51 0.73 0.81];

%% ========================================================================
%  PRINT
%% ========================================================================
print_table('MAIN. Present-value own-revenue multiplier', ...
            'numerator and denominator both cumulated to H', ...
            INS, nI, H, med_pv, lo_pv, hi_pv, N_draws, n_success);

fprintf('\n  Share of self-financing draws (total tax revenue does not fall)\n');
for j = 1:nI
    if all(isnan(sf_share(j,:))), continue; end
    fprintf('  %-22s %6.0f%% %6.0f%% %6.0f%% %6.0f%%\n', INS(j).label, 100*sf_share(j,:));
end

print_table('SECONDARY. Mechanical own-cost normalization', ...
            'no-behavior denominator; same convention as the HANK paper', ...
            INS, nI, H, med_mc, lo_mc, hi_mc, N_draws, n_success);

print_table('ROBUSTNESS. Total-revenue normalization', ...
            'denominator nets revenue induced elsewhere in the tax system', ...
            INS, nI, H, med_tt, lo_tt, hi_tt, N_draws, n_success);

fprintf('\n  Zubairy (2014) Table 2 reference values\n');
fprintf('  ---------------------------------------------------------------\n');
print_multiplier_row('Govt spending', z_med_g, z_p05_g, z_p95_g);
print_multiplier_row('Labor tax',     z_med_w, z_p05_w, z_p95_w);
print_multiplier_row('Capital tax',   z_med_k, z_p05_k, z_p95_k);

%% ========================================================================
%  WRITE CSVs
%% ========================================================================
write_mult_csv('table_multipliers_pvown.csv', INS, nI, med_pv, lo_pv, hi_pv, {});
write_mult_csv('table_multipliers_pvown_inf.csv', INS, nI, med_pi, lo_pi, hi_pi, {});
write_share_csv('table_self_financing.csv', INS, nI, sf_share);
write_mult_csv('table_multipliers_mech.csv',  INS, nI, med_mc, lo_mc, hi_mc, {});

zrows = { 'C','Govt spending', z_med_g, z_p05_g, z_p95_g ; ...
          'C','Labor tax',     z_med_w, z_p05_w, z_p95_w ; ...
          'C','Capital tax',   z_med_k, z_p05_k, z_p95_k };
write_mult_csv('table_multipliers_totalrev.csv', INS, nI, med_tt, lo_tt, hi_tt, zrows);

fprintf('\nSaved:\n');
fprintf('  table_multipliers_pvown.csv     (MAIN, finite ratios)\n');
fprintf('  table_self_financing.csv        (share of self-financing draws)\n');
fprintf('  table_multipliers_pvown_inf.csv (MAIN with the +Inf override)\n');
fprintf('  table_multipliers_mech.csv      (secondary, matches HANK paper)\n');
fprintf('  table_multipliers_totalrev.csv  (robustness, with Zubairy Panel C)\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================
function s = mk(label, panel, shock, ownvar, numsign, totsign, usetax, rho, rule)
    s = struct('label',label,'panel',panel,'shock',shock,'ownvar',ownvar, ...
               'numsign',numsign,'totsign',totsign,'usetax',usetax, ...
               'rho',rho,'rule',rule);
end

function k = idx_or_die(list, nm, what)
    k = find(strcmp(list, nm));
    if isempty(k)
        error('Could not find %s "%s".', what, nm);
    end
end

function [med, lo, hi] = summarize_block(X, clip_val)
    X = fiscal_reporting.clip(X, clip_val);
    nI = size(X,2); nH = size(X,3);
    med = NaN(nI,nH); lo = NaN(nI,nH); hi = NaN(nI,nH);
    for j = 1:nI
        for h = 1:nH
            c = X(:,j,h); c = c(~isnan(c));
            if ~isempty(c)
                stats = fiscal_reporting.percentile(c, [50 5 95]);
                med(j,h) = stats(1);
                lo(j,h)  = stats(2);
                hi(j,h)  = stats(3);
            end
        end
    end
end

function print_table(title_str, subtitle, INS, nI, H, med, lo, hi, N_draws, n_success)
    fprintf('\n=======================================================================\n');
    fprintf('%s\n', title_str);
    fprintf('  %s\n', subtitle);
    fprintf('=======================================================================\n');
    fprintf('Horizon (quarters):                   %5d  %5d  %5d  %5d\n', H);
    fprintf('(%d draws, %d successful)\n\n', N_draws, n_success);
    panels = {'A','B'};
    for pp = 1:2
        if strcmp(panels{pp},'A')
            fprintf('Panel A. Zubairy-original instruments\n');
        else
            fprintf('\nPanel B. Extension instruments\n');
        end
        for j = 1:nI
            if strcmp(INS(j).panel, panels{pp})
                print_multiplier_row(INS(j).label, med(j,:), lo(j,:), hi(j,:));
            end
        end
    end
end

function write_mult_csv(fname, INS, nI, med, lo, hi, extra_rows)
    rows = {};
    for j = 1:nI
        rows(end+1,:) = [{INS(j).panel, INS(j).label, 'median'}, ...
                         arrayfun(@(v) sprintf('%.2f',v), med(j,:), 'UniformOutput', false)]; %#ok<AGROW>
        rows(end+1,:) = [{INS(j).panel, INS(j).label, '[5,95]'}, ...
                         arrayfun(@(a,b) sprintf('[%.2f,%.2f]',a,b), lo(j,:), hi(j,:), ...
                         'UniformOutput', false)]; %#ok<AGROW>
    end
    for r = 1:size(extra_rows,1)
        rows(end+1,:) = [{extra_rows{r,1}, extra_rows{r,2}, 'median'}, ...
                         arrayfun(@(v) sprintf('%.2f',v), extra_rows{r,3}, 'UniformOutput', false)]; %#ok<AGROW>
        rows(end+1,:) = [{extra_rows{r,1}, extra_rows{r,2}, '[5,95]'}, ...
                         arrayfun(@(a,b) sprintf('[%.2f,%.2f]',a,b), extra_rows{r,4}, extra_rows{r,5}, ...
                         'UniformOutput', false)]; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'Panel','Multiplier','Statistic','H1','H4','H12','H20'});
    writetable(T, fname, 'Delimiter', ',', 'QuoteStrings', true);
end

function write_share_csv(fname, INS, nI, share)
    rows = {};
    for j = 1:nI
        if all(isnan(share(j,:))), continue; end
        rows(end+1,:) = [{INS(j).panel, INS(j).label, 'share'}, ...
                         arrayfun(@(v) sprintf('%.3f',v), share(j,:), 'UniformOutput', false)]; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'Panel','Multiplier','Statistic','H1','H4','H12','H20'});
    writetable(T, fname, 'Delimiter', ',', 'QuoteStrings', true);
end

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
            if ~exist(mh_file, 'file'), continue; end
            S = load(mh_file, 'x2');
            if ~isfield(S, 'x2') || isempty(S.x2), continue; end
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

function out = safe_ratio(num, den, net_cost)
    out = fiscal_reporting.multiplier(num, den, net_cost);
end
