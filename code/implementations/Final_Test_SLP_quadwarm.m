
rng('shuffle');
num_runs = 50;
timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

save_starts = false;            % Save generated starts only when use_saved_starts is false.
use_saved_starts = true;        % Use saved starts if true; generate new starts if false.


saved_starts_file = 'fould3_starting_point';   % Name of the saved starting-point file.



problem_library = [
    make_problem("fould2",  "foulds2.nl",     @generate_foulds2_start,  "low", -1100)
    make_problem("fould3",  "foulds3.nl",     @generate_fould3_start,  "low", -8)
    make_problem("fould4",  "foulds4.nl",      @generate_fould4_start,  "low",-8)
    make_problem("fould5",  "foulds5.nl",     @generate_fould5_start,  "low", -8)
    make_problem("bental4", "bental4new.nl",   @generate_bental4_start, "low", -450)
    make_problem("bental5", "bental5_new.nl",  @generate_bental5_start, "low", -3500)
    make_problem("adhya1",  "adhya1.nl",     @generate_adhya1_start,  "low", -549.8030502)
    make_problem("adhya2",  "adhya2.nl",     @generate_adhya2_start,  "low", -549.8030502)
    make_problem("adhya3",  "adhya3.nl",     @generate_adhya3_start,  "low", -561.0446875)
    make_problem("adhya4",  "adhya4.nl",     @generate_adhya4_start,  "low", -877.6457399)
    make_problem("rt2",     "rt2.nl",        @generate_rt2_start,     "low", -4391.825894)
];


selected_problem_names = ["fould3"];   % Select the problem name from problem_library above.
if any(selected_problem_names == "all")
    problems = problem_library;
else
    problems = problem_library(ismember([problem_library.name], selected_problem_names));
end
if isempty(problems)
    error('No configured problem matched selected_problem_names.');
end

for problem_id = 1:numel(problems)
    cfg = problems(problem_id);
    prob = NonlinearProblem(char(cfg.nlfile));
    n = get_nvar(prob);
    m = get_ncon(prob);
    bl = get_bl(prob);
    bu = get_bu(prob);
    cl = get_cl(prob);
    cu = get_cu(prob);

    if use_saved_starts
        S = load(saved_starts_file, 'starts');
        starts = S.starts;

        if size(starts, 1) ~= n
            error('Saved starts dimension does not match this problem.');
        end

        num_runs_this = size(starts, 2);
    else
        starts = NaN(n, num_runs);
        num_runs_this = num_runs;
    end

    batch_results = repmat(make_empty_result(n), num_runs_this, 1);


    fprintf('\n=== SLP quadprog warm | %s | %s | mode=%s | runs=%d ===\n', ...
      char(cfg.name), char(cfg.nlfile), char(cfg.start_mode), num_runs_this);
    fprintf('%5s %9s %14s %12s %10s %8s %8s %10s %18s %24s\n', ...
      'run', 'success', 'final_f', 'final_h', 'iters', 'rest', 'maxviol', 'time_s', 'status', 'reason');

    for run_id = 1:num_runs_this
        run_timer = tic;
        try
            if use_saved_starts
                x_start = starts(:, run_id);
            else
                x_start = make_start(cfg);
                starts(:, run_id) = x_start;
            end

            [result, ~] = solve_slp_quadwarm_once(prob, n, m, bl, bu, cl, cu, x_start, cfg);
        catch ME
            result = make_runtime_error_result(ME, n);
        end
        result.solve_time = toc(run_timer);


        batch_results(run_id) = result;

        fprintf('%5d %9d %14.6g %12.3e %10d %8d %8d %10.4f %18s %24s\n', ...
            run_id, result.success, result.final_f, result.final_h, ...
            result.iterations, result.num_restoration, result.max_viol_idx, ...
            result.solve_time, char(result.status), char(result.failure_reason));
    end

    summary = summarize_batch(batch_results, cfg);
    disp('=== batch summary ===')
    disp(summary)



    method_name = "SLP_warm";
    result_file = char(method_name + "_" + cfg.name + "_" + timestamp + ".mat");

    save(result_file, ...
        'batch_results', 'summary', 'starts', ...
        'saved_starts_file', '-v7');

    fprintf('Saved results to %s\n', result_file);


    if save_starts && ~use_saved_starts
        starts_file = char("SLP_quadwarm_starts_" + cfg.name + "_" + timestamp + ".mat");
        save(starts_file, "starts");
        fprintf('Saved starts to %s\n', starts_file);
    end
end


function cfg = make_problem(name, nlfile, start_fun, start_mode, target_f)
cfg = struct();
cfg.name = string(name);
cfg.nlfile = string(nlfile);
cfg.start_fun = start_fun;
cfg.start_mode = string(start_mode);
cfg.target_f = target_f;
cfg.target_tol = 1e-3;
cfg.feas_tol = 1e-5;
cfg.max_iter = 200;
end


function x = make_start(cfg)
try
    x = cfg.start_fun(cfg.start_mode);
catch
    x = cfg.start_fun();
end
x = x(:);
end


function summary = summarize_batch(batch_results, cfg)
success_vec = [batch_results.success];
f_vec = [batch_results.final_f];
h_vec = [batch_results.final_h];
time_vec = [batch_results.solve_time];
global_success = success_vec & abs(f_vec - cfg.target_f) < cfg.target_tol & h_vec < cfg.feas_tol;
feasible_local = ~global_success & h_vec < cfg.feas_tol;
near_feasible = ~global_success & h_vec >= cfg.feas_tol & h_vec < 1e-4;
failed = h_vec >= 1e-4;
finite_f = f_vec(isfinite(f_vec));

if isempty(finite_f)
    best_objective = NaN;
    worst_objective = NaN;
    best_gap_to_target = NaN;
else
    best_objective = min(finite_f);
    worst_objective = max(finite_f);
    best_gap_to_target = min(abs(finite_f - cfg.target_f));
end

summary = struct( ...
    'num_runs', numel(batch_results), ...
    'global_success_count', sum(global_success), ...
    'feasible_local_count', sum(feasible_local), ...
    'near_feasible_count', sum(near_feasible), ...
    'failed_count', sum(failed), ...
    'best_objective', best_objective, ...
    'worst_objective', worst_objective, ...
    'best_gap_to_target', best_gap_to_target, ...
    'max_final_h', max(h_vec),...
    'mean_time', mean(time_vec, 'omitnan'), ...
    'total_time', sum(time_vec, 'omitnan'), ...
    'max_time', max(time_vec));
end


function result = make_empty_result(n)
result = struct();
result.success = false;
result.final_f = NaN;
result.final_h = Inf;
result.max_viol = Inf;
result.max_viol_idx = NaN;
result.iterations = 0;
result.x0 = NaN(n,1);
result.x_final = NaN(n,1);
result.f0 = NaN;
result.h0 = Inf;
result.num_restoration = 0;
result.status = "";
result.failure_reason = "";
result.solve_time = NaN;
end


function result = make_runtime_error_result(ME, n)
result = make_empty_result(n);
result.status = "runtime_error";
result.failure_reason = string(ME.message);
end


function [result, log] = solve_slp_quadwarm_once(prob, n, m, bl, bu, cl, cu, x, cfg)

% Record the initial point.
x0 = x;
f0 = eval_obj(prob, x);
c0 = eval_cons(prob, 1:m, x);
[h0, viol0] = constraint_violation(c0, cl, cu);
[max_viol0, max_viol0_idx] = max(viol0);

log = struct();
log.x0 = x0;
log.f0 = f0;
log.h0 = h0;
log.max_viol0 = max_viol0;
log.max_viol0_idx = max_viol0_idx;
log.iter = struct([]);
log.restoration = struct([]);
log.restoration_soc = struct([]);
log.main_soc = struct([]);
log.final = struct();
log.failure_reason = "";

last_main_rejected = false;
rho = 1;
rho_min = 1e-8;
max_iter = cfg.max_iter;

filter_f = [];
filter_h = [];
c_init = eval_cons(prob, 1:m, x);
[h_init, ~] = constraint_violation(c_init, cl, cu);
filter_upper_bound = max(0.5, 1.25 * h_init);
first_after_restoration = false;
qp_options = optimoptions('quadprog', ...
    'Display','none', ...
    'Algorithm','active-set', ...
    'ConstraintTolerance',1e-6, ...
    'OptimalityTolerance',1e-6, ...
    'MaxIterations',5000);

small_progress_count = 0;
small_progress_limit = 10;
ws_main = [];

for k = 1:max_iter
    rho_before = rho;
    x_current = x;
    f = eval_obj(prob, x_current);
    g = eval_obj_grad(prob, x_current);
    cval = eval_cons(prob, 1:m, x_current);
    J = eval_jac_val(prob, 1:m, x_current);
    J = J';
    [Acanon, ccanon] = canonical_linearization(cval, J, cl, cu);

    % Record the main outer iteration.
    [h_current, ~] = constraint_violation(cval, cl, cu);
    log.iter(k).k = k;
    log.iter(k).f = f;
    log.iter(k).h_current = h_current;
    log.iter(k).rho = rho;
    log.iter(k).norm_g_inf = norm(g, inf);
    log.iter(k).entered_restoration = false;
    log.iter(k).main_qp_exitflag = NaN;
    log.iter(k).norm_d_inf = NaN;
    log.iter(k).accepted = false;
    log.iter(k).f_trial = NaN;
    log.iter(k).h_trial = NaN;
    log.iter(k).pred_obj_change = NaN;
    log.iter(k).actual_obj_change = NaN;

    if h_current < 1e-4 && rho < 1e-6
        if ~isfield(log, 'rho_reset_count')
            log.rho_reset_count = 0;
        end
        if log.rho_reset_count < 3
            rho = 1e-3;
            log.rho_reset_count = log.rho_reset_count + 1;
        else
            log.failure_reason = "near_feasible_stagnation";
            break;
        end
    end

    % Detect long-term infeasibility stagnation to avoid empty iterations up to max_iter.
    if k > 10
        recent_h = [log.iter(k-9:k).h_current];
        if max(recent_h) - min(recent_h) < 1e-5 && h_current > 1e-4
            log.failure_reason = "stagnation_infeasible";
            break;
        end
    end

    Aqp = Acanon;
    bqp = -ccanon;
    lb_d = max(bl - x_current, -rho * ones(n,1));
    ub_d = min(bu - x_current,  rho * ones(n,1));
    restoration_phase = false;

    Hqp = zeros(n);
    if isempty(ws_main)
        ws_main = optimwarmstart(zeros(n,1), qp_options);
    end

    [ws_main, ~, exitflag, ~, lambda] = quadprog_warm_silent( ...
        Hqp, g(:), Aqp, bqp, [], [], lb_d, ub_d, ws_main);
    d = ws_main.X;

    % Record the main LP/QP result.
    log.iter(k).main_qp_exitflag = exitflag;
    if exist('d', 'var') && ~isempty(d)
        log.iter(k).norm_d_inf = norm(d, inf);
    end

    if exitflag <= 0
        restoration_phase = true;
        log.iter(k).entered_restoration = true;
        max_rest_iter = 20;
        restoration_tol = 1e-8;

        x_rest = x_current;
        restoration_success = false;
        rest_filter_hJ = [];
        rest_filter_hJperp = [];
        ccanon_start = canonical_value(cval, cl, cu);
        rest_upper_bound = max(sum(max(ccanon_start, 0)), restoration_tol);
        prev_idxJ = [];
        best_iter_x = x_rest;
        best_iter_h_l1 = rest_upper_bound;

        % Record restoration-phase entry.
        log.iter(k).rest_start_l1 = sum(max(ccanon_start, 0));
        log.iter(k).rest_upper_bound = rest_upper_bound;
        log.iter(k).restoration_success = false;

        for rest_iter = 1:max_rest_iter
            c_rest = eval_cons(prob, 1:m, x_rest);
            J_rest = eval_jac_val(prob, 1:m, x_rest);
            J_rest = J_rest';
            [Acanon_rest, ccanon_rest] = canonical_linearization(c_rest, J_rest, cl, cu);

            lb_d_rest = max(bl - x_rest, -rho * ones(n,1));
            ub_d_rest = min(bu - x_rest,  rho * ones(n,1));

            [~, s_elastic, lp1_feasible, idxJ, idxJperp, elasticflag, elasticfval, elasticoutput] = ...
                solve_elastic_lp_partition( ...
                Acanon_rest, ccanon_rest, lb_d_rest, ub_d_rest, ...
                n, restoration_tol, qp_options);

            rest_log_id = numel(log.restoration) + 1;
            [h_rest_now, ~] = constraint_violation(c_rest, cl, cu);
            h_l1_now = sum(max(ccanon_rest, 0));

            log.restoration(rest_log_id).h_l1_now = h_l1_now;
            log.restoration(rest_log_id).outer_iter = k;
            log.restoration(rest_log_id).rest_iter = rest_iter;
            log.restoration(rest_log_id).h_rest = h_rest_now;
            log.restoration(rest_log_id).elasticflag = elasticflag;
            log.restoration(rest_log_id).elasticfval = elasticfval;
            if isempty(s_elastic)
                log.restoration(rest_log_id).total_slack = NaN;
            else
                log.restoration(rest_log_id).total_slack = sum(s_elastic);
            end
            log.restoration(rest_log_id).lp1_feasible = lp1_feasible;
            log.restoration(rest_log_id).num_idxJ = sum(idxJ);
            log.restoration(rest_log_id).restflag = NaN;
            log.restoration(rest_log_id).norm_drest_inf = NaN;
            log.restoration(rest_log_id).h_l1_trial = NaN;
            log.restoration(rest_log_id).pred_l1_trial = NaN;
            log.restoration(rest_log_id).true_vs_pred = NaN;
            log.restoration(rest_log_id).true_decrease = false;
            log.restoration(rest_log_id).pred_decrease = false;
            log.restoration(rest_log_id).h_l1_ratio = NaN;
            log.restoration(rest_log_id).h_rest_ratio = NaN;
            log.restoration(rest_log_id).accepted_rest = false;

            if elasticflag <= 0
                break;
            end

            if lp1_feasible
                g_rest = eval_obj_grad(prob, x_rest);
                Hmain = zeros(n);
                dmain0 = zeros(n,1);
                [dmain, ~, exitflag, ~, lambda] = quadprog_silent( ...
                    Hmain, g_rest(:), Acanon_rest, -ccanon_rest, [], [], ...
                    lb_d_rest, ub_d_rest, dmain0, qp_options);

                if exitflag > 0
                    x_main_trial = x_rest + dmain;
                    f_main_trial = eval_obj(prob, x_main_trial);
                    c_main_trial = eval_cons(prob, 1:m, x_main_trial);
                    [h_main_trial, ~] = constraint_violation(c_main_trial, cl, cu);

                    obj_worsen_limit = 0.05 * max(1, abs(f));
                    switch_objective_not_destroyed = f_main_trial <= f + obj_worsen_limit;

                    if h_current > 1e-1
                        switch_objective_guard_ok = true;
                    else
                        switch_objective_guard_ok = switch_objective_not_destroyed;
                    end

                    main_switch_ok = h_main_trial <= h_rest_now && ...
                        filter_accept(filter_f, filter_h, f_main_trial, h_main_trial) && ...
                        switch_objective_guard_ok;

                    if main_switch_ok
                        d = x_main_trial - x_current;
                        restoration_phase = false;
                        break;
                    end
                end
            end

            if h_rest_now < max(1e-6, 0.01 * h_current)
                restoration_success = true;
                break;
            end

            if ~any(idxJ)
                idxJ = ccanon_rest > restoration_tol;
                idxJperp = ~idxJ;
            end
            if isempty(prev_idxJ)
                J_changed = false;
            else
                J_changed = any(idxJ ~= prev_idxJ);
            end
            prev_idxJ = idxJ;
            nJperp = sum(idxJperp);

            [drest, restflag, ~, restfval, restoutput] = solve_phase1_lp( ...
                Acanon_rest, ccanon_rest, idxJ, idxJperp, ...
                lb_d_rest, ub_d_rest, n, qp_options);

            log.restoration(rest_log_id).restflag = restflag;
            log.restoration(rest_log_id).restfval = restfval;
            log.restoration(rest_log_id).rest_message = string(restoutput.message);
            if isfield(restoutput, 'iterations')
                log.restoration(rest_log_id).rest_iterations = restoutput.iterations;
            else
                log.restoration(rest_log_id).rest_iterations = NaN;
            end
            if isfield(restoutput, 'funcCount')
                log.restoration(rest_log_id).rest_funcCount = restoutput.funcCount;
            else
                log.restoration(rest_log_id).rest_funcCount = NaN;
            end

            if exist('drest', 'var') && ~isempty(drest)
                rest_lp_resid = Acanon_rest(idxJperp,:) * drest + ccanon_rest(idxJperp);
                if isempty(rest_lp_resid)
                    log.restoration(rest_log_id).rest_constrviolation = 0;
                else
                    log.restoration(rest_log_id).rest_constrviolation = max(max(rest_lp_resid, 0));
                end
            else
                log.restoration(rest_log_id).rest_constrviolation = NaN;
            end

            if isfield(restoutput, 'firstorderopt')
                log.restoration(rest_log_id).rest_firstorderopt = restoutput.firstorderopt;
            else
                log.restoration(rest_log_id).rest_firstorderopt = NaN;
            end
            if exist('drest', 'var') && ~isempty(drest)
                log.restoration(rest_log_id).norm_drest_inf = norm(drest, inf);
            end

            if restflag <= 0
                rest_constrviol = Inf;
                if exist('drest', 'var') && ~isempty(drest)
                    rest_lp_resid = Acanon_rest(idxJperp,:) * drest + ccanon_rest(idxJperp);
                    if isempty(rest_lp_resid)
                        rest_constrviol = 0;
                    else
                        rest_constrviol = max(max(rest_lp_resid, 0));
                    end
                end

                allow_near_feasible_rest_lp = ...
                    restflag == -2 && ...
                    isscalar(rest_constrviol) && ...
                    isfinite(rest_constrviol) && ...
                    rest_constrviol <= 1e-5 && ...
                    exist('drest', 'var') && ...
                    ~isempty(drest);

                if ~allow_near_feasible_rest_lp
                    break;
                end
            end

            x_rest_trial = x_rest + drest;
            c_rest_trial = eval_cons(prob, 1:m, x_rest_trial);
            [hJ_trial, hJperp_trial, h_l1_trial] = ...
                phaseI_violation_values(c_rest_trial, cl, cu, idxJ, idxJperp);
            [h_rest_trial, ~] = constraint_violation(c_rest_trial, cl, cu);

            pred_l1_trial = sum(max(ccanon_rest + Acanon_rest * drest, 0));
            true_vs_pred = h_l1_trial / max(pred_l1_trial, 1e-12);
            true_decrease = h_l1_trial < h_l1_now;
            pred_decrease = pred_l1_trial < h_l1_now;

            best_iter_x = x_rest_trial;
            best_iter_h = h_rest_trial;
            best_iter_h_l1 = h_l1_trial;

            filter_ok = filter_accept(rest_filter_hJ, rest_filter_hJperp, hJ_trial, hJperp_trial);
            l1_decrease = h_l1_trial < h_l1_now;
            max_not_worse = h_rest_trial <= h_rest_now;
            accepted_rest = filter_ok && l1_decrease && max_not_worse;

            log.restoration(rest_log_id).hJ_trial = hJ_trial;
            log.restoration(rest_log_id).hJperp_trial = hJperp_trial;
            log.restoration(rest_log_id).h_l1_trial = h_l1_trial;
            log.restoration(rest_log_id).h_rest_trial = h_rest_trial;
            log.restoration(rest_log_id).pred_l1_trial = pred_l1_trial;
            log.restoration(rest_log_id).true_vs_pred = true_vs_pred;
            log.restoration(rest_log_id).true_decrease = true_decrease;
            log.restoration(rest_log_id).pred_decrease = pred_decrease;
            log.restoration(rest_log_id).h_l1_ratio = h_l1_trial / max(h_l1_now, 1e-12);
            log.restoration(rest_log_id).h_rest_ratio = h_rest_trial / max(h_rest_now, 1e-12);
            log.restoration(rest_log_id).accepted_rest = accepted_rest;

            stop_restoration = false;
            if accepted_rest
                x_rest = x_rest_trial;
                restoration_success = true;
                [rest_filter_hJ, rest_filter_hJperp] = filter_update( ...
                    rest_filter_hJ, rest_filter_hJperp, hJ_trial, hJperp_trial);

                if h_l1_trial < rest_upper_bound
                    rest_upper_bound = max(h_l1_trial, rest_upper_bound / 10);
                end

                if ~last_main_rejected && norm(drest,inf) >= 0.8 * rho
                    rho = min(2*rho, 50);
                end
                if h_rest_trial < 1e-6
                    stop_restoration = true;
                end
            else
                rho = max(0.5 * rho, rho_min);
            end
            if stop_restoration
                break;
            end
        end

        if restoration_phase
            if restoration_success
                d = x_rest - x_current;
                log.iter(k).restoration_success = true;
            else
                log.failure_reason = "restoration_failed";
                break;
            end
        end
    end

    x_trial = x_current + d;
    f_trial = eval_obj(prob, x_trial);
    c_trial = eval_cons(prob, 1:m, x_trial);
    [h_current, ~] = constraint_violation(cval, cl, cu);
    [h_trial, ~] = constraint_violation(c_trial, cl, cu);

    % Record the main trial point.
    log.iter(k).f_trial = f_trial;
    log.iter(k).h_trial = h_trial;
    log.iter(k).pred_obj_change = g' * d;
    log.iter(k).actual_obj_change = f_trial - f;
    if exist('d', 'var') && ~isempty(d)
        log.iter(k).norm_d_inf = norm(d, inf);
    end

    if restoration_phase
        accepted = h_trial < h_current;
    else
        accepted = h_trial <= filter_upper_bound && filter_accept(filter_f,filter_h, f_trial, h_trial);
        if first_after_restoration && h_trial > h_current
            accepted = false;
        end
    end

    if accepted
        log.iter(k).accepted = true;
        x = x_trial;
        last_main_rejected = false;

        if ~restoration_phase
            [filter_f, filter_h] = filter_update(filter_f, filter_h, f_trial, h_trial);
        end
        if ~restoration_phase && norm(d,inf) >= 0.8 * rho_before
            rho = min(1.5*rho, 10);
        end
        if log.iter(k).entered_restoration
            first_after_restoration = true;
        else
            first_after_restoration = false;
        end
    else
        log.iter(k).accepted = false;
        rho = max(0.5 * rho, rho_min);
        first_after_restoration = false;
        last_main_rejected = true;
    end

    if h_trial < 1e-5 && abs(f_trial - f) < 1e-6
        small_progress_count = small_progress_count + 1;
    else
        small_progress_count = 0;
    end

    if small_progress_count >= small_progress_limit
        log.failure_reason = "small_progress";
        break;
    end

    if norm(d,inf) < 1e-6 && h_trial < 1e-5
        log.failure_reason = "small_step";
        break;
    end
end

% Record final results.
c_final = eval_cons(prob, 1:m, x);
[h_final, viol_final] = constraint_violation(c_final, cl, cu);
f_final = eval_obj(prob, x);
[max_viol_final, max_viol_final_idx] = max(viol_final);
log.final.x = x;
log.final.f = f_final;
log.final.h = h_final;
log.final.max_viol = max_viol_final;
log.final.max_viol_idx = max_viol_final_idx;
log.final.iterations = k;
log.final.success = h_final < cfg.feas_tol && abs(f_final - cfg.target_f) < cfg.target_tol;

if ~log.final.success && strlength(log.failure_reason) == 0
    if h_final < cfg.feas_tol
        log.failure_reason = "feasible_local_solution";
    elseif h_final < 1e-4
        log.failure_reason = "near_feasible_not_global";
    elseif k >= max_iter
        log.failure_reason = "max_iter_reached";
    else
        log.failure_reason = "not_converged";
    end
end

result = make_empty_result(n);
result.success = log.final.success;
result.final_f = f_final;
result.final_h = h_final;
result.max_viol = max_viol_final;
result.max_viol_idx = max_viol_final_idx;
result.iterations = k;
result.x0 = x0;
result.x_final = x;
result.f0 = f0;
result.h0 = h0;
result.num_restoration = numel(log.restoration);
result.status = classify_run_status(log.final.success, h_final, log.failure_reason, cfg);
result.failure_reason = log.failure_reason;
end


function [h, viol] = constraint_violation(c, cl, cu)
viol_lower = max(cl - c, 0);
viol_upper = max(c - cu, 0);
viol_lower(~isfinite(cl)) = 0;
viol_upper(~isfinite(cu)) = 0;
viol = max(viol_lower, viol_upper);
h = max(viol);
end


function [Acanon, ccanon] = canonical_linearization(c, J, cl, cu)
idxU = isfinite(cu);
idxL = isfinite(cl);
Acanon = [J(idxU,:); -J(idxL,:)];
ccanon = [c(idxU) - cu(idxU); cl(idxL) - c(idxL)];
end


function ccanon = canonical_value(c, cl, cu)
idxU = isfinite(cu);
idxL = isfinite(cl);
ccanon = [c(idxU) - cu(idxU); cl(idxL) - c(idxL)];
end


function accepted = filter_accept(filter_a, filter_b, a_new, b_new)
accepted = true;
for i = 1:length(filter_a)
    if a_new >= filter_a(i) && b_new >= filter_b(i)
        accepted = false;
        return
    end
end
end


function [filter_a, filter_b] = filter_update(filter_a, filter_b, a_new, b_new)
keep = true(length(filter_a),1);
for i = 1:length(filter_a)
    if filter_a(i) >= a_new && filter_b(i) >= b_new
        keep(i) = false;
    end
end
filter_a = filter_a(keep);
filter_b = filter_b(keep);
filter_a = [filter_a; a_new];
filter_b = [filter_b; b_new];
end


function [dlp, s_elastic, lp1_feasible, idxJ, idxJperp, elasticflag, elasticfval, elasticoutput] = ...
    solve_elastic_lp_partition(Acanon, ccanon, lb_d, ub_d, n, restoration_tol, qp_options)

ncanon = length(ccanon);
felastic = [zeros(n,1); ones(ncanon,1)];
Aelastic = [Acanon, -eye(ncanon)];
belastic = -ccanon;
lbelastic = [lb_d; zeros(ncanon,1)];
ubelastic = [ub_d; Inf(ncanon,1)];

Hz = zeros(n + ncanon);
z0 = [zeros(n,1); max(ccanon, 0)];
z0 = min(max(z0, lbelastic), ubelastic);

[zelastic, elasticfval, elasticflag, elasticoutput] = quadprog_silent( ...
    Hz, felastic, Aelastic, belastic, [], [], ...
    lbelastic, ubelastic, z0, qp_options);

if elasticflag <= 0
    dlp = [];
    s_elastic = [];
    lp1_feasible = false;
    idxJ = [];
    idxJperp = [];
    return
end

dlp = zelastic(1:n);
s_elastic = zelastic(n+1:end);
total_slack = sum(s_elastic);
lp1_feasible = total_slack <= restoration_tol;
J_tol = max(restoration_tol, 1e-5);
idxJ = s_elastic > J_tol;
idxJperp = ~idxJ;
end


function [drest, restflag, restfun, restfval, restoutput] = solve_phase1_lp( ...
    Acanon_rest, ccanon_rest, idxJ, idxJperp, ...
    lb_d_rest, ub_d_rest, n, qp_options)

if any(idxJ)
    frest = sum(Acanon_rest(idxJ,:), 1)';
else
    frest = zeros(n,1);
end

Arest = [];
brest = [];

if any(idxJperp)
    Arest = Acanon_rest(idxJperp,:);
    brest = -ccanon_rest(idxJperp);
end

restfun = @(drest) frest' * drest;

Hrest = zeros(n);
drest0 = zeros(n,1);

[drest, restfval, restflag, restoutput] = quadprog_silent( ...
    Hrest, frest, Arest, brest, [], [], ...
    lb_d_rest, ub_d_rest, drest0, qp_options);
end

function [x, fval, exitflag, output, lambda] = quadprog_silent( ...
    H, f, A, b, Aeq, beq, lb, ub, x0, options)
warning_state = warning('off', 'all');
cleanup = onCleanup(@() warning(warning_state));
[x, fval, exitflag, output, lambda] = quadprog( ...
    H, f, A, b, Aeq, beq, lb, ub, x0, options);
end


function [ws, fval, exitflag, output, lambda] = quadprog_warm_silent( ...
    H, f, A, b, Aeq, beq, lb, ub, ws)
warning_state = warning('off', 'all');
cleanup = onCleanup(@() warning(warning_state));
[ws, fval, exitflag, output, lambda] = quadprog( ...
    H, f, A, b, Aeq, beq, lb, ub, ws);
end


function [hJ, hJperp, h_l1] = phaseI_violation_values(c, cl, cu, idxJ, idxJperp)
ccanon = canonical_value(c, cl, cu);
hJ = sum(max(ccanon(idxJ), 0));
hJperp = sum(max(ccanon(idxJperp), 0));
h_l1 = sum(max(ccanon, 0));
end

function status = classify_run_status(success, h_final, failure_reason, cfg)
if success
    status = "global";
elseif h_final < cfg.feas_tol
    status = "local";
elseif h_final < 1e-4
    status = "near_feasible";
elseif strlength(failure_reason) > 0
    status = failure_reason;
else
    status = "not_converged";
end
end
