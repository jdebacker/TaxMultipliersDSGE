function run_uk(what)
% run_uk  Estimation driver for tax_dsge_uk (2026-10-24).
%
%   run_uk('check')     parse, steady state, Blanchard-Kahn / rank check, and
%                       start-file names vs estimated objects (about a minute)
%   run_uk('estimate')  check, then newrat mode search from start_uk.mat
%                       (run 2: UK run-1 best draw and chain covariance, ITC dropped), then
%                       MCMC 4 x 1M
%   run_uk('post')      reload chains (STAGE 5), compute_multipliers,
%                       post_estimation
%
% Log: tax_dsge_uk\run_uk_log.txt

if nargin < 1, what = 'check'; end
D   = fileparts(mfilename('fullpath'));
setappdata(0, 'UK_LOG', fullfile(D, 'run_uk_log.txt'));
say('==================== run_uk : %s ====================', what);

switch what
    case 'check'
        preflight(D);
    case 'estimate'
        if ~preflight(D), say('PREFLIGHT FAILED. Nothing estimated.'); return; end
        ok1 = stage1(D);
        stage2(D, ok1);
    case 'post'
        if stage_run(D, 5, 'reload chains (no sampling)')
            for scr = {'compute_multipliers', 'post_estimation'}
                t = tic; say('  %s', scr{1});
                try
                    cd(D); evalin('base', scr{1});
                    say('    done in %.1f min', toc(t)/60);
                catch ME
                    say('    FAIL: %s', ME.message);
                end
            end
        end
    otherwise
        error('Unknown option "%s". Use check, estimate or post.', what);
end
cd(D); say('run_uk %s DONE.', what);
end

% =========================================================================
function ok = preflight(D)
global M_ estim_params_ %#ok<GVMIS>
ok = false;
say('PREFLIGHT %s', D);
try
    cd(D);
    if ~exist('start_uk.mat', 'file'), error('start_uk.mat missing'); end
    dynare tax_dsge_uk noclearall          % STAGE undefined -> steady + check only
    names = [cellstr(M_.exo_names(estim_params_.var_exo(:,1))); ...
             cellstr(M_.param_names(estim_params_.param_vals(:,1)))];
    S = load('start_uk.mat');
    sn = strtrim(cellstr(S.parameter_names)); sn = sn(:);
    if numel(sn) ~= numel(names) || ~all(strcmp(sn, names(:)))
        error('start file names do not match estimated objects (%d vs %d)', numel(sn), numel(names));
    end
    pv = @(n) M_.params(strcmp(cellstr(M_.param_names), n));
    say('  OK  %d estimated objects match start file; tauw %.3f tauk %.3f tauc %.3f tauic %.3f g_y %.3f b_y %.2f', ...
        numel(names), pv('tauw_ss'), pv('tauk_ss'), pv('tauc_ss'), pv('tauic_ss'), pv('g_y'), pv('b_y'));
    ok = true;
catch ME
    say('  FAIL %s', ME.message);
end
end

% =========================================================================
function ok = stage1(D)
ok = false; t = tic; t0 = now;
say('STAGE 1 newrat mode search');
try
    cd(D);
    copyfile('start_uk.mat', 'start_mode.mat', 'f');
    dynare tax_dsge_uk noclearall -DSTAGE=1
    mf = fullfile('tax_dsge_uk', 'Output', 'tax_dsge_uk_mode.mat');
    d = dir(mf);
    if isempty(d) || d.datenum < t0, error('newrat did not write a new %s', mf); end
    copyfile(mf, 'mode_newrat.mat', 'f');
    O = load('mode_newrat.mat');
    say('  done in %.0f min. UK log posterior at mode %.2f', toc(t)/60, -O.fval);
    ok = true;
catch ME
    say('  FAIL after %.0f min: %s', toc(t)/60, ME.message);
end
end

% =========================================================================
function stage2(D, ok1)
t = tic;
say('STAGE 2 MCMC 4 x 1M');
try
    cd(D);
    S = load('start_mode.mat');
    if ok1
        O = load('mode_newrat.mat');
        on = strtrim(cellstr(O.parameter_names)); sn = strtrim(cellstr(S.parameter_names));
        if ~isequal(on(:), sn(:)), error('mode and start names differ'); end
        xparam1 = O.xparam1; fval = O.fval;
        say('  mode/fval from newrat (%.2f); proposal = chain covariance from start_uk.mat', -fval);
    else
        xparam1 = S.xparam1; fval = S.fval;
        say('  WARNING newrat failed: starting at the US best draw; fval is from the US model, so C4 is NOT valid.');
    end
    hh = S.hh; parameter_names = strtrim(cellstr(S.parameter_names)); %#ok<NASGU>
    save('tax_dsge_uk_mode_cal.mat', 'xparam1', 'hh', 'fval', 'parameter_names');
    dynare tax_dsge_uk noclearall -DSTAGE=2
    say('  MCMC done in %.1f h', toc(t)/3600);
catch ME
    say('  FAIL after %.1f h: %s', toc(t)/3600, ME.message);
end
end

% =========================================================================
function ok = stage_run(D, stage, label)
ok = false; t = tic;
say('STAGE %d %s', stage, label);
try
    cd(D);
    eval(sprintf('dynare tax_dsge_uk noclearall -DSTAGE=%d', stage));
    say('  done in %.1f min', toc(t)/60);
    ok = true;
catch ME
    say('  FAIL after %.1f min: %s', toc(t)/60, ME.message);
end
end

% =========================================================================
function say(varargin)
msg = sprintf('[%s] %s', datestr(now, 'yyyy-mm-dd HH:MM:SS'), sprintf(varargin{:}));
fprintf('%s\n', msg);
LOG = getappdata(0, 'UK_LOG');
if ~isempty(LOG)
    fid = fopen(LOG, 'a'); if fid > 0, fprintf(fid, '%s\r\n', msg); fclose(fid); end
end
end
