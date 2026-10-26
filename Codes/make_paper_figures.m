function make_paper_figures(what)
% make_paper_figures  Figures and posterior tables for the draft (2026-10-26).
%
%   make_paper_figures          compute from the chains of the five final
%                               models, then plot
%   make_paper_figures('plot')  replot from paper_figures\fig_data.mat only
%
% Chains of tax_dsge, tax_dsge_ar1, tax_dsge_kgoods, tax_dsge_kgoods_ar1
% and tax_dsge_uk are read READ-ONLY. Dynare is run without estimation on a
% copy of each .mod in paper_figures\_work\<model>\, so nothing inside the
% model folders is written.
%
% Draws: the same 1000 posterior draws as compute_multipliers.m (rng 2014,
% mh_drop 0.5), PHI_ZERO = true.
% Multipliers: present-value own-revenue (main convention), H = 1..40.
% IRFs: one-standard-deviation expansionary shock (tax cut; increase in
% spending, ITC or expensing). Output, consumption, investment and hours in
% percent deviation from steady state; inflation and the nominal rate in
% annualized percentage points. Median and 5-95 percentile band.
%
% Output, in Codes\paper_figures\:
%   irf_eps_g/w/k/c/ic/etau.png  main model (tax_dsge)
%   mult_main.png                main model multipliers by horizon, 90% band
%   mult_models.png              four US models, posterior medians
%   mult_uk.png                  UK vs US AR(1)
%   fig_data.mat, mult_paths.csv, posterior_<model>.csv, log.txt
%
% The log compares the H1/H4/H12/H20 medians with each folder's
% table_multipliers_pvown.csv (should agree to rounding).

if nargin < 1, what = 'all'; end
C   = fileparts(mfilename('fullpath'));          % ...\Codes
OUT = fullfile(C, 'paper_figures');
if ~exist(OUT, 'dir'), mkdir(OUT); end
setappdata(0, 'PF_LOG', fullfile(OUT, 'log.txt'));
say('==================== make_paper_figures : %s ====================', what);

% folder, fname, data file, label
MODELS = { 'tax_dsge',            'tax_dsge_est', 'US_Data_Matlab_1961_2013.mat', 'Full feedback'
           'tax_dsge_ar1',        'tax_dsge_est', 'US_Data_Matlab_1961_2013.mat', 'AR(1)'
           'tax_dsge_kgoods',     'tax_dsge_est', 'US_Data_Matlab_1961_2013.mat', 'Full feedback + friction'
           'tax_dsge_kgoods_ar1', 'tax_dsge_est', 'US_Data_Matlab_1961_2013.mat', 'AR(1) + friction'
           'tax_dsge_uk',         'tax_dsge_uk',  'UK_Data_Matlab_tax.mat',       'UK, AR(1)' };

% label, shock, own revenue/outlay, sign on PV(dY), IRF direction, label in table csv
INS = { 'Government spending', 'eps_g',    'g',        +1, +1, 'Govt spending'
        'Labor tax',           'eps_w',    'rev_w',    -1, -1, 'Labor tax'
        'Capital tax',         'eps_k',    'rev_k',    -1, -1, 'Capital tax'
        'Consumption tax',     'eps_c',    'rev_c',    -1, -1, 'Consumption tax'
        'ITC',                 'eps_ic',   'spend_ic', +1, +1, 'ITC'
        'Expensing',           'eps_etau', 'rev_k',    -1, +1, 'Expensing rate' };

% variable, label, scale type ('pct' = 100*dx/xss, 'app' = 400*dx)
VARS = { 'y',  'Output',        'pct'
         'c',  'Consumption',   'pct'
         'i',  'Investment',    'pct'
         'h',  'Hours',         'pct'
         'pi', 'Inflation',     'app'
         'R',  'Interest rate', 'app' };
T = 40;

if strcmp(what, 'plot')
    L = load(fullfile(OUT, 'fig_data.mat'));
    R = L.R;
else
    R = [];
    for m = 1:size(MODELS, 1)
        r = run_model(C, OUT, MODELS(m,:), INS, VARS, T);
        if isempty(r), say('STOPPED: %s failed. Nothing plotted.', MODELS{m,1}); cd(C); return; end
        R = [R r]; %#ok<AGROW>
    end
    save(fullfile(OUT, 'fig_data.mat'), 'R', 'MODELS', 'INS', 'VARS', 'T');
    write_paths(OUT, R, INS, T);
end
plot_all(OUT, R, INS, VARS, T);
cd(C);
say('make_paper_figures DONE. Output in %s', OUT);
end

% =========================================================================
function r = run_model(C, OUT, spec, INS, VARS, T)
global M_ oo_ options_ estim_params_ %#ok<GVMIS>
r = [];
folder = spec{1}; fname = spec{2}; t0 = tic;
say('%s', folder);
try
    src = fullfile(C, folder);
    W   = fullfile(OUT, '_work', folder);
    if ~exist(W, 'dir'), mkdir(W); end
    copyfile(fullfile(src, [fname '.mod']), W, 'f');
    copyfile(fullfile(src, spec{3}), W, 'f');
    cd(W);
    eval(sprintf('dynare %s noclearall', fname));   % no -DSTAGE: steady state and check only
    % Stay in W until the simulations are done: the steady-state and dynamic
    % files Dynare generated (package +<fname>) are found from this folder.
    M0 = M_; oo0 = oo_; opt0 = options_; ep = estim_params_;
    % With no estimation or stoch_simul command, qz_criterium stays empty and
    % order stays at Dynare's default; resol then fails on every draw.
    if ~isfield(opt0, 'qz_criterium') || isempty(opt0.qz_criterium), opt0.qz_criterium = 1 + 1e-6; end
    opt0.order = 1;   % Dynare's default is 2 when no estimation command sets it
    % Fields that the estimation command normally fills in (same fix as the
    % epsK script): object counts used by set_all_parameters, and the state-
    % space ordering used by resol.
    ep.nvx = size(ep.var_exo,1); ep.nvn = size(ep.var_endo,1); ep.ncx = size(ep.corrx,1);
    ep.ncn = size(ep.corrn,1);   ep.np  = size(ep.param_vals,1);
    if ~isfield(oo0, 'dr') || ~isfield(oo0.dr, 'kstate')
        if ~isfield(oo0, 'dr'), oo0.dr = struct(); end
        try
            oo0.dr = set_state_space(oo0.dr, M0, opt0);
        catch
            oo0.dr = set_state_space(oo0.dr, M0);
        end
    end
    if ~strcmp(M0.fname, fname), error('loaded model is %s, expected %s', M0.fname, fname); end
    pn = cellstr(M0.param_names); vn = cellstr(M0.endo_names); sn = cellstr(M0.exo_names);
    has_kg = any(strcmp(pn, 'epsK'));
    if has_kg ~= contains(folder, 'kgoods'), error('friction parameter epsK present = %d in %s', has_kg, folder); end

    % ----- draws: identical selection to compute_multipliers.m -----
    MH = fullfile(src, fname, 'metropolis');
    S  = load(fullfile(MH, [fname '_mh_history_0.mat']), 'record');
    rng(2014);
    post = read_draws(MH, fname, S.record.Nblck, S.record.LastFileNumber, 0.5);
    nest = size(ep.var_exo,1) + size(ep.var_endo,1) + size(ep.corrx,1) + size(ep.corrn,1) + size(ep.param_vals,1);
    if size(post, 2) ~= nest, error('draw length %d vs %d estimated objects', size(post,2), nest); end
    n_post = size(post, 1);
    say('  %d retained draws, %d estimated objects', n_post, nest);

    % ----- posterior table (all retained draws) -----
    if ~isempty(ep.var_endo) || ~isempty(ep.corrx) || ~isempty(ep.corrn)
        error('posterior table assumes only shock std devs and parameters are estimated');
    end
    names = [strcat('stderr_', cellstr(M0.exo_names(ep.var_exo(:,1)))); cellstr(M0.param_names(ep.param_vals(:,1)))];
    pri   = [ep.var_exo(:, 5:7); ep.param_vals(:, 5:7)];
    shp   = {'beta', 'gamma', 'normal', 'invgamma', 'uniform', 'invgamma2'};
    pshape = cell(nest, 1);
    for k = 1:nest, pshape{k} = shp{pri(k,1)}; end
    Tpost = table(names, pshape, pri(:,2), pri(:,3), mean(post)', median(post)', prctile(post, 5)', prctile(post, 95)', ...
        'VariableNames', {'name', 'prior', 'prior_mean', 'prior_sd', 'post_mean', 'post_median', 'p05', 'p95'});
    writetable(Tpost, fullfile(OUT, sprintf('posterior_%s.csv', folder)));

    try
        M_t = set_all_parameters(mean(post)', ep, M0);
        [~, info_t] = resol(0, M_t, opt0, oo0);
        if info_t(1) ~= 0, error('resol info = %d', info_t(1)); end
    catch MEt
        error('test solve at the posterior mean failed: %s (%s line %d)', MEt.message, MEt.stack(1).name, MEt.stack(1).line);
    end
    say('  test solve at the posterior mean OK');
    N = min(1000, n_post);
    sel = post(randperm(n_post, N), :);
    clear post

    % ----- indices -----
    nI = size(INS, 1); nV = size(VARS, 1);
    have = false(nI, 1); si = zeros(nI, 1); vo = zeros(nI, 1);
    for j = 1:nI
        a = find(strcmp(sn, INS{j,2})); b = find(strcmp(vn, INS{j,3}));
        if ~isempty(a) && ~isempty(b), have(j) = true; si(j) = a; vo(j) = b; end
    end
    vv = zeros(nV, 1);
    for v = 1:nV
        k = find(strcmp(vn, VARS{v,1})); if isempty(k), error('variable %s not found', VARS{v,1}); end
        vv(v) = k;
    end
    iy = find(strcmp(vn, 'y')); iR = find(strcmp(vn, 'R')); ip_phi = find(strcmp(pn, 'phi'));

    % ----- simulate -----
    MULT = NaN(N, nI, T); IRF = NaN(N, nI, nV, T); ok = false(N, 1); first_err = '';
    for d = 1:N
        try
            M_d = set_all_parameters(sel(d,:)', ep, M0);
            if ~isempty(ip_phi), M_d.params(ip_phi) = 0; end
            oo_d = oo0;
            [dr, info, M_d, oo_d] = resol(0, M_d, opt0, oo_d);
            if info(1) ~= 0
                if isempty(first_err), first_err = sprintf('resol info = %d', info(1)); end
                continue
            end
            ys = oo_d.dr.ys; disc = ys(iR) .^ (-(0:T-1)');
            for j = find(have)'
                ex = zeros(T, M_d.exo_nbr); ex(1, si(j)) = 1;
                yp = simult_(M_d, opt0, ys, dr, ex, 1);
                if size(yp, 2) == T + 1, irf = yp(:, 2:end) - ys; else, irf = yp - ys; end
                den = cumsum(disc .* irf(vo(j),:)');
                mlt = INS{j,4} * cumsum(disc .* irf(iy,:)') ./ den;
                mlt(abs(den) < 1e-10 | abs(mlt) > 50) = NaN;
                MULT(d, j, :) = mlt;
                sd = sqrt(M_d.Sigma_e(si(j), si(j)));
                for v = 1:nV
                    if strcmp(VARS{v,3}, 'pct'), x = 100 * irf(vv(v),:) / ys(vv(v)); else, x = 400 * irf(vv(v),:); end
                    IRF(d, j, v, :) = INS{j,5} * sd * x;
                end
            end
            ok(d) = true;
        catch ME2
            if isempty(first_err), first_err = sprintf('%s (%s line %d)', ME2.message, ME2.stack(1).name, ME2.stack(1).line); end
        end
        if d == 20 && ~any(ok(1:20)), error('first 20 draws failed; first error: %s', first_err); end
    end
    cd(C);
    say('  solved %d of %d draws', sum(ok), N);
    if ~isempty(first_err), say('  first failure: %s', first_err); end
    if sum(ok) < 0.9*N, error('only %d of %d draws solved', sum(ok), N); end

    r.folder = folder; r.label = spec{4}; r.have = have; r.nsolved = sum(ok);
    r.mult_med = squeeze(median(MULT(ok,:,:), 1, 'omitnan'));
    r.mult_lo  = squeeze(prctile(MULT(ok,:,:), 5, 1));
    r.mult_hi  = squeeze(prctile(MULT(ok,:,:), 95, 1));
    r.irf_med  = squeeze(median(IRF(ok,:,:,:), 1, 'omitnan'));
    r.irf_lo   = squeeze(prctile(IRF(ok,:,:,:), 5, 1));
    r.irf_hi   = squeeze(prctile(IRF(ok,:,:,:), 95, 1));

    % ----- check against the folder's table -----
    try
        Tt = readtable(fullfile(src, 'table_multipliers_pvown.csv'), 'TextType', 'string', 'Delimiter', ',');
        Hc = [1 4 12 20]; worst = 0;
        for j = find(have)'
            row = Tt(Tt.Multiplier == INS{j,6} & Tt.Statistic == "median", :);
            if isempty(row), continue; end
            tab = [num(row.H1) num(row.H4) num(row.H12) num(row.H20)];
            dd = abs(tab - r.mult_med(j, Hc)); if any(isnan(dd)), worst = NaN; break; end
            worst = max(worst, max(dd));
        end
        say('  max |median - table| at H1/H4/H12/H20 = %.3f (table is rounded to 0.01)', worst);
    catch ME
        say('  table check skipped: %s', ME.message);
    end
    for j = find(have)'
        say('  %-20s H1 %5.2f  H12 %5.2f  H20 %5.2f  H40 %5.2f', INS{j,1}, r.mult_med(j,1), r.mult_med(j,12), r.mult_med(j,20), r.mult_med(j,40));
    end
    say('  done in %.1f min', toc(t0)/60);
catch ME
    cd(C);
    say('  FAIL: %s', ME.message);
    r = [];
end
end

% =========================================================================
function write_paths(OUT, R, INS, T)
rows = {};
for m = 1:numel(R)
    for j = find(R(m).have)'
        for h = 1:T
            rows(end+1, :) = {R(m).folder, INS{j,1}, h, R(m).mult_med(j,h), R(m).mult_lo(j,h), R(m).mult_hi(j,h)}; %#ok<AGROW>
        end
    end
end
writetable(cell2table(rows, 'VariableNames', {'model', 'instrument', 'H', 'median', 'p05', 'p95'}), fullfile(OUT, 'mult_paths.csv'));
end

% =========================================================================
function plot_all(OUT, R, INS, VARS, T)
t = (1:T)';
fol = {R.folder};
main = R(strcmp(fol, 'tax_dsge'));
COL = [0.00 0.30 0.65; 0.85 0.33 0.10; 0.00 0.30 0.65; 0.85 0.33 0.10];
STY = {'-', '-', '--', '--'};
band = [0.75 0.85 0.95];

% ----- IRFs, main model -----
for j = find(main.have)'
    f = newfig(900, 560);
    for v = 1:size(VARS, 1)
        subplot(2, 3, v);
        shade(t, squeeze(main.irf_lo(j,v,:)), squeeze(main.irf_hi(j,v,:)), band);
        plot(t, squeeze(main.irf_med(j,v,:)), '-', 'Color', COL(1,:), 'LineWidth', 1.6);
        finish(VARS{v,2}, T, v > 3, ifelse(strcmp(VARS{v,3}, 'pct'), '% dev.', 'ann. pp'));
    end
    savefig_png(f, fullfile(OUT, sprintf('irf_%s.png', INS{j,2})));
end

% ----- multipliers by horizon, main model -----
f = newfig(900, 560);
k = 0;
for j = find(main.have)'
    k = k + 1; subplot(2, 3, k);
    shade(t, main.mult_lo(j,:)', main.mult_hi(j,:)', band);
    plot(t, main.mult_med(j,:)', '-', 'Color', COL(1,:), 'LineWidth', 1.6);
    finish(INS{j,1}, T, k > 3, '');
end
savefig_png(f, fullfile(OUT, 'mult_main.png'));

% ----- four US models -----
US = {'tax_dsge', 'tax_dsge_ar1', 'tax_dsge_kgoods', 'tax_dsge_kgoods_ar1'};
f = newfig(900, 560); h = gobjects(1, 4); ax3 = [];
for j = 1:size(INS, 1)
    ax = subplot(2, 3, j); hold on;
    for m = 1:4
        r = R(strcmp(fol, US{m}));
        hh = plot(t, r.mult_med(j,:)', STY{m}, 'Color', COL(m,:), 'LineWidth', 1.5);
        if j == 3, h(m) = hh; end
    end
    if j == 3, ax3 = ax; end
    finish(INS{j,1}, T, j > 3, '');
end
% legend in the capital tax panel, whose lower part is empty
legend(ax3, h, cellfun(@(u) R(strcmp(fol, u)).label, US, 'UniformOutput', false), 'Location', 'southeast', 'FontSize', 7);
savefig_png(f, fullfile(OUT, 'mult_models.png'));

% ----- UK vs US AR(1) -----
uk = R(strcmp(fol, 'tax_dsge_uk')); us = R(strcmp(fol, 'tax_dsge_ar1'));
f = newfig(900, 560); k = 0;
for j = find(uk.have)'
    k = k + 1; subplot(2, 3, k);
    shade(t, uk.mult_lo(j,:)', uk.mult_hi(j,:)', band);
    h1 = plot(t, uk.mult_med(j,:)', '-', 'Color', COL(1,:), 'LineWidth', 1.6);
    h2 = plot(t, us.mult_med(j,:)', '--', 'Color', [0 0 0], 'LineWidth', 1.3);
    finish(INS{j,1}, T, k > 2, '');
    if k == 1, legend([h1 h2], {'UK', 'US AR(1)'}, 'Location', 'best', 'FontSize', 8); end
end
savefig_png(f, fullfile(OUT, 'mult_uk.png'));
end

% =========================================================================
function f = newfig(w, h)
f = figure('Position', [100 100 w h], 'Visible', 'off', 'Color', 'white', 'PaperPositionMode', 'auto');
end

function shade(t, lo, hi, c)
fill([t; flipud(t)], [lo; flipud(hi)], c, 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on;
end

function finish(ttl, T, bottom, ylab)
plot([1 T], [0 0], ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6);
title(ttl, 'FontSize', 10, 'FontWeight', 'bold');
if bottom, xlabel('Quarters', 'FontSize', 9); end
if ~isempty(ylab), ylabel(ylab, 'FontSize', 8); end
xlim([1 T]); grid on; set(gca, 'FontSize', 8, 'Box', 'on'); hold off;
end

function savefig_png(f, fn)
print(f, fn, '-dpng', '-r300'); close(f); say('  saved %s', fn);
end

function x = num(v)
if isnumeric(v), x = double(v(1)); else, x = str2double(string(v(1))); end
end

function o = ifelse(c, a, b)
if c, o = a; else, o = b; end
end

% =========================================================================
function P = read_draws(mhpath, root, nb, nf, drop)
% Same as post_estimation.m / compute_multipliers.m.
P = [];
for b = 1:nb
    Cb = [];
    for f = 1:nf
        fn = fullfile(mhpath, sprintf('%s_mh%d_blck%d.mat', root, f, b));
        if ~exist(fn, 'file'), continue; end
        S = load(fn, 'x2');
        if isfield(S, 'x2') && ~isempty(S.x2), Cb = [Cb; S.x2]; end %#ok<AGROW>
    end
    if isempty(Cb), warning('No draws for block %d.', b); continue; end
    P = [P; Cb(ceil(drop*size(Cb,1))+1:end, :)]; %#ok<AGROW>
end
end

function say(varargin)
msg = sprintf('[%s] %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), sprintf(varargin{:}));
fprintf('%s\n', msg);
LOG = getappdata(0, 'PF_LOG');
if ~isempty(LOG)
    fid = fopen(LOG, 'a'); if fid > 0, fprintf(fid, '%s\r\n', msg); fclose(fid); end
end
end
