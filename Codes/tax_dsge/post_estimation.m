% post_estimation.m
% =========================================================================
% Post-estimation robustness on the converged chains (2026-10-23).
%
% Same 1000 posterior draws as compute_multipliers.m (rng 2014, same burn-in),
% PHI_ZERO = true, present-value own-revenue multipliers (MAIN convention).
% For every draw the model is re-solved under each setting below with all
% estimated parameters held at that draw, and the multiplier is compared
% draw by draw with the baseline.
%
%   baseline            the calibration as estimated (sig_u 2.5, bg 0.70,
%                       epsK 4 where the capital-goods block exists)
%   sig_u = 2.0, 3.0    utilization curvature (drives the capital tax)
%   bg    = 0.55, 0.85  government deep habit (drives spending)
%   epsK  = 6, 9, 14    capital-goods supply elasticity, House-Shapiro range
%   epsK  = 1e6         no friction (block collapses to the no-friction model)
%
% These are conditional re-solves, not re-estimations: the posterior is the
% one obtained at sig_u = 2.5, bg = 0.70, epsK = 4.
%
% Run (from this folder):
%   dynare tax_dsge_est noclearall -DSTAGE=5
%   post_estimation
%
% Writes: post_pvown.csv   median and [5,95] by setting, instrument, horizon
%         post_paired.csv  paired % change vs baseline, median and [5,95]
%         post_raw.mat
% =========================================================================

if ~exist('M_','var') || ~exist('oo_','var') || ~exist('options_','var') || ~exist('estim_params_','var')
    error('Run Dynare first: dynare tax_dsge_est noclearall -DSTAGE=5');
end
if ~strcmp(M_.fname,'tax_dsge_est')
    error('Loaded workspace is "%s", expected "tax_dsge_est".', M_.fname);
end

N_draws  = 1000;
H        = [1 4 12 20];
T_irf    = max(H);
clip_val = 50;
PHI_ZERO = true;

pn = cellstr(M_.param_names); vn = cellstr(M_.endo_names); sn = cellstr(M_.exo_names);
fi = @(list,nm) find(strcmp(list,nm));
has_kg = ~isempty(fi(pn,'epsK'));

% ----- settings: {label, parameter, value}; first row is the baseline -----
SET = { 'baseline', '', NaN
        'sig_u=2.0', 'sig_u', 2.0
        'sig_u=3.0', 'sig_u', 3.0
        'bg=0.55',   'bg',    0.55
        'bg=0.85',   'bg',    0.85 };
if has_kg
    SET = [SET; { 'epsK=6','epsK',6; 'epsK=9','epsK',9; 'epsK=14','epsK',14; 'epsK=1e6 (no friction)','epsK',1e6 }];
end
nS = size(SET,1);
ipS = zeros(nS,1);
for e = 2:nS
    k = fi(pn, SET{e,2}); if isempty(k), error('parameter %s not found', SET{e,2}); end
    ipS(e) = k;
end
ip_phi = fi(pn,'phi');

% ----- instruments (interest income tax is exogenous: not included) -----
INS = { 'Govt spending',   'eps_g',   'g',        +1
        'Labor tax',       'eps_w',   'rev_w',    -1
        'Capital tax',     'eps_k',   'rev_k',    -1
        'Transfers',       'eps_tr',  'tr',       +1
        'Consumption tax', 'eps_c',   'rev_c',    -1
        'ITC',             'eps_ic',  'spend_ic', +1
        'Dividend tax',    'eps_td',  'rev_d',    -1
        'Expensing rate',  'eps_etau','rev_k',    -1 };
nI = size(INS,1); nH = numel(H);
si = zeros(nI,1); vi = zeros(nI,1);
for j = 1:nI
    a = fi(sn, INS{j,2}); b = fi(vn, INS{j,3});
    if isempty(a) || isempty(b), error('instrument %s not found', INS{j,1}); end
    si(j) = a; vi(j) = b;
end
ns = cell2mat(INS(:,4)); iy = fi(vn,'y'); iR = fi(vn,'R');

% ----- draws: identical selection to compute_multipliers.m -----
if isfield(options_,'mh_drop') && ~isempty(options_.mh_drop), mh_drop = options_.mh_drop; else, mh_drop = 0.333; end
rng(2014);
MH = fullfile(M_.dname, 'metropolis');
load(fullfile(MH, [M_.fname '_mh_history_0.mat']), 'record');
post = read_draws(MH, M_.fname, record.Nblck, record.LastFileNumber, mh_drop);
n_post = size(post,1);
nest = size(estim_params_.var_exo,1) + size(estim_params_.var_endo,1) + size(estim_params_.corrx,1) ...
     + size(estim_params_.corrn,1) + size(estim_params_.param_vals,1);
if size(post,2) ~= nest, error('draw length %d vs %d estimated objects', size(post,2), nest); end
if n_post < N_draws, N_draws = n_post; end
sel = post(randperm(n_post, N_draws), :);
fprintf('\npost_estimation: %d draws, %d settings, %d instruments, PHI_ZERO = %d\n', N_draws, nS, nI, PHI_ZERO);

M_base = M_; oo_base = oo_; opt_base = options_;
NUM = NaN(nS, N_draws, nI, nH);  DEN = NUM;  ok = false(nS, N_draws);
t0 = tic;
for d = 1:N_draws
    M_d0 = set_all_parameters(sel(d,:)', estim_params_, M_base);
    if PHI_ZERO, M_d0.params(ip_phi) = 0; end
    for e = 1:nS
        try
            M_d = M_d0;
            if e > 1, M_d.params(ipS(e)) = SET{e,3}; end
            oo_d = oo_base;
            [dr, info, M_d, oo_d] = resol(0, M_d, opt_base, oo_d);
            if info(1) ~= 0, continue; end
            ys = oo_d.dr.ys; disc = ys(iR) .^ (-(0:T_irf-1)');
            for j = 1:nI
                ex = zeros(T_irf, M_d.exo_nbr); ex(1, si(j)) = 1;
                yp = simult_(M_d, opt_base, ys, dr, ex, 1);
                if size(yp,2) == T_irf + 1, irf = yp(:,2:end) - ys; else, irf = yp - ys; end
                pv_y = ns(j) * cumsum(disc .* irf(iy,:)');
                pv_o = cumsum(disc .* irf(vi(j),:)');
                NUM(e,d,j,:) = pv_y(H);  DEN(e,d,j,:) = pv_o(H);
            end
            ok(e,d) = true;
        catch
        end
    end
    if mod(d,200) == 0, fprintf('  %4d / %d  solved per setting: %s  (%.0f s)\n', d, N_draws, mat2str(sum(ok,2)'), toc(t0)); end
end
M_ = M_base; oo_ = oo_base; options_ = opt_base;

PV = NUM ./ DEN; PV(abs(DEN) < 1e-10) = NaN; PV(abs(PV) > clip_val) = NaN;
both = all(ok,1);
fprintf('\nSolved: %s of %d draws per setting; paired sample %d\n', mat2str(sum(ok,2)'), N_draws, sum(both));

% ----- table 1: levels -----
rows = {};
for e = 1:nS
    for j = 1:nI
        labs = {'median','p05','p95'};
        for st = 1:3
            v = NaN(1,nH);
            for h = 1:nH
                x = squeeze(PV(e,ok(e,:),j,h)); x = x(~isnan(x));
                if isempty(x), continue; end
                if st == 1, v(h) = median(x); elseif st == 2, v(h) = prctile(x,5); else, v(h) = prctile(x,95); end
            end
            rows(end+1,:) = [{SET{e,1}, INS{j,1}, labs{st}}, num2cell(v)]; %#ok<AGROW>
        end
    end
end
writetable(cell2table(rows,'VariableNames',{'Setting','Instrument','Statistic','H1','H4','H12','H20'}), 'post_pvown.csv');

% ----- table 2: paired % change vs baseline -----
rows = {};
for e = 2:nS
    for j = 1:nI
        for h = 1:nH
            a = squeeze(PV(e,both,j,h)); b = squeeze(PV(1,both,j,h));
            p = 100*(a./b - 1); p = p(isfinite(p));
            if isempty(p), m = NaN; lo = NaN; hi = NaN; else, m = median(p); lo = prctile(p,5); hi = prctile(p,95); end
            rows(end+1,:) = {SET{e,1}, INS{j,1}, H(h), median(b,'omitnan'), median(a,'omitnan'), m, lo, hi}; %#ok<AGROW>
        end
    end
end
writetable(cell2table(rows,'VariableNames',{'Setting','Instrument','H','m_baseline','m_setting','pct_chg_median','pct_p05','pct_p95'}), 'post_paired.csv');

% ----- screen: H12 and H20 medians -----
for hh = [3 4]
    fprintf('\nH%d pv-own medians (paired sample)\n%-17s', H(hh), '');
    for e = 1:nS, fprintf(' %9s', strtok(SET{e,1},' ')); end; fprintf('\n');
    for j = 1:nI
        fprintf('%-17s', INS{j,1});
        for e = 1:nS, fprintf(' %9.3f', median(squeeze(PV(e,both,j,hh)),'omitnan')); end
        fprintf('\n');
    end
end
save('post_raw.mat','SET','INS','H','NUM','DEN','ok','sel','PHI_ZERO');
fprintf('\nSaved post_pvown.csv, post_paired.csv, post_raw.mat  (%.1f min)\n', toc(t0)/60);

% -------------------------------------------------------------------------
function P = read_draws(mhpath, root, nb, nf, drop)
    P = [];
    for b = 1:nb
        C = [];
        for f = 1:nf
            fn = fullfile(mhpath, sprintf('%s_mh%d_blck%d.mat', root, f, b));
            if ~exist(fn,'file'), continue; end
            S = load(fn, 'x2');
            if isfield(S,'x2') && ~isempty(S.x2), C = [C; S.x2]; end %#ok<AGROW>
        end
        if isempty(C), warning('No draws for block %d.', b); continue; end
        P = [P; C(ceil(drop*size(C,1))+1:end, :)]; %#ok<AGROW>
    end
end
