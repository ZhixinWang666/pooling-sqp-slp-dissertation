


rng('shuffle');
num_runs = 50;


timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));

save_starts = false;            % Save generated starts only when use_saved_starts is false.
use_saved_starts =  true;       % Use saved starts if true; generate new starts if false.


saved_starts_file = 'fould5_starting_point';   % Name of the saved starting-point file.

problem_library = [
    make_problem("fould2",  "foulds2.nl",     @generate_foulds2_start,  "low", -1100)
    make_problem("fould3",  "foulds3.nl",     @generate_fould3_start,  "low", -8)
    make_problem("fould4",  "foulds4.nl",      @generate_fould4_start,  "low",-8)
    make_problem("fould5",  "foulds5.nl",     @generate_fould5_start,  "low", -8)
    make_problem("bental4", "bental4new.nl",     @generate_bental4_start, "low", -450)
    make_problem("bental5", "bental5_new.nl",     @generate_bental5_start, "low", -3500)
    make_problem("adhya1",  "adhya1.nl",     @generate_adhya1_start,  "low", -549.8030502)
    make_problem("adhya2",  "adhya2.nl",     @generate_adhya2_start,  "low", -549.8030502)
    make_problem("adhya3",  "adhya3.nl",     @generate_adhya3_start,  "low", -561.0446875)
    make_problem("adhya4",  "adhya4.nl",     @generate_adhya4_start,  "low", -877.6457399)
    make_problem("rt2",     "rt2.nl",        @generate_rt2_start,     "low",    -4391.825894)
];


selected_problem_names = ["fould5"];   % Select the problem name from problem_library above.
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
   

    fprintf('\n=== %s | %s | mode=%s | runs=%d ===\n', ...
      char(cfg.name), char(cfg.nlfile), char(cfg.start_mode), num_runs_this);
    fprintf('%5s %9s %14s %12s %10s %8s %8s %8s %10s %18s %24s\n', ...
        'run', 'success', 'final_f', 'final_h', 'iters', 'rest', 'soc', 'maxviol', 'time_s', 'status','reason');


    for run_id = 1:num_runs_this
        run_timer = tic;
        try
            if use_saved_starts
                x_start = starts(:, run_id);
            else
                x_start = cfg.start_fun(cfg.start_mode);
                x_start = x_start(:);
                starts(:, run_id) = x_start;
            end

            [result, ~] = solve_pool_once(prob, n, m, bl, bu, cl, cu, x_start, cfg);
        catch ME
            result = make_runtime_error_result(ME, n);
            
        end
        result.solve_time = toc(run_timer);

        batch_results(run_id) = result;
        
            
        

        fprintf('%5d %9d %14.6g %12.3e %10d %8d %8d %8d %10.4f %18s %24s\n', ...
            run_id, result.success, result.final_f, result.final_h, ...
            result.iterations, result.num_restoration, result.num_main_soc, ...
            result.max_viol_idx, result.solve_time, char(result.status), char(result.failure_reason));
    end

    summary = summarize_batch(batch_results, cfg);
    disp('=== batch summary ===')
    disp(summary)

    method_name = "SQP_old";
    result_file = char(method_name + "_" + cfg.name + "_" + timestamp + ".mat");

    save(result_file, ...
        'batch_results', 'summary', 'starts', ...
        'saved_starts_file', '-v7');

    fprintf('Saved results to %s\n', result_file);


    
    
       
       
    
    if save_starts && ~use_saved_starts
        starts_file = char("SQP_starts_" + cfg.name + "_" + timestamp + ".mat");
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
cfg.filter_upper_base = 0.5;
cfg.rho0 = 1;
end


function summary = summarize_batch(batch_results, cfg)
success_vec = [batch_results.success];
f_vec = [batch_results.final_f];
h_vec = [batch_results.final_h];
time_vec = [batch_results.solve_time];
global_success = success_vec & abs(f_vec - cfg.target_f) < cfg.target_tol & h_vec <cfg.feas_tol;
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
    'max_final_h', max(h_vec), ...
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
result.num_restoration_soc = 0;
result.num_main_soc = 0;
result.status = "";
result.failure_reason = "";
result.failed_rest_message = "";
result.failed_restflag = NaN;
result.failed_rest_constrviolation = NaN;
result.failed_rest_firstorderopt = NaN;
result.solve_time = NaN;
end


function result = make_runtime_error_result(ME, n)
result = make_empty_result(n);
result.status = "runtime_error";
result.failure_reason = string(ME.message);
end


function [result, log] = solve_pool_once(prob, n, m, bl, bu, cl, cu, x, cfg)

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
% End initial-point record.

lambda_current = zeros(m,1);
rho = cfg.rho0;
rho_min = 1e-8;
max_iter = cfg.max_iter;
rest_hessian_mode = "positive_lagrangian";
filter_f = [];
filter_h = [];
c_init = eval_cons(prob, 1:m, x);
[h_init, ~] = constraint_violation(c_init, cl, cu); 

filter_upper_bound = max(cfg.filter_upper_base, 1.25 * h_init);
first_after_restoration = false;
last_main_rejected = false;
small_progress_count = 0;
small_progress_limit = 10;
options = optimoptions('fmincon', ...
    'Display','none', ...
    'Algorithm','sqp', ...
    'SpecifyObjectiveGradient',true, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-6, ...
    'ConstraintTolerance',1e-6);
for k = 1:max_iter
    rho_before = rho;
    x_current = x;
    lambda_before = lambda_current;
    f = eval_obj(prob, x_current);
    g = eval_obj_grad(prob, x_current);
    cval = eval_cons(prob, 1:m, x_current);
    J = eval_jac_val(prob, 1:m, x_current);
    J = J';
    H = eval_hesslag(prob, lambda_before, x_current);
    H = 0.5 * (H + H');
    [Acanon, ccanon] = canonical_linearization(cval, J, cl, cu);

    % Record the main outer iteration.
    [h_current, ~] = constraint_violation(cval, cl, cu);
    minEigH = min(eig(H));
    maxCanon = max(max(ccanon, 0));
    log.iter(k).k = k;
    log.iter(k).f = f;
    log.iter(k).h_current = h_current;
    log.iter(k).rho = rho;
    log.iter(k).norm_g_inf = norm(g, inf);
    log.iter(k).minEigH = minEigH;
    log.iter(k).maxCanon = maxCanon;
    log.iter(k).entered_restoration = false;
    log.iter(k).main_qp_exitflag = NaN;
    log.iter(k).norm_d_inf = NaN;
    log.iter(k).accepted = false;
    log.iter(k).f_trial = NaN;
    log.iter(k).h_trial = NaN;
    % End main outer-iteration record.
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
    qpfun = @(d) quad_objective(d, H, g);
    d0 = zeros(n,1);
    [d, ~, exitflag, ~, lambda] = fmincon( ...
        qpfun, d0, Aqp, bqp, [], [], lb_d, ub_d, [], options);

    % Record the main QP result.
    log.iter(k).main_qp_exitflag = exitflag;
    if exist('d', 'var') && ~isempty(d)
        log.iter(k).norm_d_inf = norm(d, inf);
    end
    % End main QP result record.

    if exitflag <= 0
        restoration_phase = true;
        log.iter(k).entered_restoration = true;
        max_rest_iter = 20;
        restoration_tol = 1e-8;
        lp_options = optimoptions('linprog','Display','none');
        x_rest = x_current;
        restoration_success = false;
        rest_filter_hJ = [];
        rest_filter_hJperp = [];
        ccanon_start = canonical_value(cval, cl, cu);
        rest_upper_bound = max(sum(max(ccanon_start, 0)), restoration_tol);
        prev_idxJ = [];
        best_iter_x = x_rest;
        best_iter_h = Inf;
        best_iter_h_l1 = rest_upper_bound;

        % Record restoration-phase entry.
        log.iter(k).rest_start_l1 = sum(max(ccanon_start, 0));
        log.iter(k).rest_upper_bound = rest_upper_bound;
        log.iter(k).restoration_success = false;
        % End restoration-phase entry record.

        for rest_iter = 1:max_rest_iter
            c_rest = eval_cons(prob, 1:m, x_rest);
            J_rest = eval_jac_val(prob, 1:m, x_rest);
            J_rest = J_rest';
            [Acanon_rest, ccanon_rest] = canonical_linearization(c_rest, J_rest, cl, cu);
            lb_d_rest = max(bl - x_rest, -rho * ones(n,1));
            ub_d_rest = min(bu - x_rest,  rho * ones(n,1));

            [dlp, s_elastic, lp1_feasible, idxJ, idxJperp, elasticflag] = ...
                solve_elastic_lp_partition( ...
                Acanon_rest, ccanon_rest, lb_d_rest, ub_d_rest, ...
                n, restoration_tol, lp_options);

            % Record the elastic LP result.
            rest_log_id = numel(log.restoration) + 1;
            [h_rest_now, ~] = constraint_violation(c_rest, cl, cu);
            h_l1_now = sum(max(ccanon_rest, 0));
            log.restoration(rest_log_id).h_l1_now = h_l1_now;
            log.restoration(rest_log_id).outer_iter = k;
            log.restoration(rest_log_id).rest_iter = rest_iter;
            log.restoration(rest_log_id).h_rest = h_rest_now;
            log.restoration(rest_log_id).elasticflag = elasticflag;
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
            log.restoration(rest_log_id).accepted_rest = false;
            % End elastic LP result record.

            if elasticflag <= 0
                log.failure_reason = "elastic_lp_failed";
                break;
            end

            if lp1_feasible
                f_rest = eval_obj(prob, x_rest);
                g_rest = eval_obj_grad(prob, x_rest);
                H_rest_main = eval_hesslag(prob, lambda_before, x_rest);
                H_rest_main = 0.5 * (H_rest_main + H_rest_main');
                qpfun_rest = @(dmain) quad_objective(dmain, H_rest_main, g_rest);
                [dmain, ~, exitflag, ~, lambda] = fmincon( ...
                    qpfun_rest, dlp, Acanon_rest, -ccanon_rest, [], [], ...
                    lb_d_rest, ub_d_rest, [], options);
                if exitflag > 0
                    % Check that the main step is acceptable before switching back from restoration.
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

            [drest, restflag, restfun, restfval, restoutput] = solve_phase1_qp( ...
                prob, x_rest, lambda_before, Acanon_rest, ccanon_rest, ...
                idxJ, idxJperp, lb_d_rest, ub_d_rest, n, ...
                rest_hessian_mode, options);

            % Record the restoration QP result.
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
            if isfield(restoutput, 'constrviolation')
                log.restoration(rest_log_id).rest_constrviolation = restoutput.constrviolation;
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
            % End restoration QP result record.

            if restflag <= 0
                rest_constrviol = Inf;
                if isfield(restoutput, 'constrviolation')
                    rest_constrviol = restoutput.constrviolation;
                end

                allow_near_feasible_rest_qp = ...
                    restflag == -2 && ...
                    rest_constrviol <= 1e-5 && ...
                    exist('drest', 'var') && ...
                    ~isempty(drest);

                if ~allow_near_feasible_rest_qp
                    warning('Restoration QP failed at outer iteration %d, restoration iteration %d', ...
                        k, rest_iter);
                    break;
                end
            end
            x_rest_trial = x_rest + drest;
            c_rest_trial = eval_cons(prob, 1:m, x_rest_trial);
            [hJ_trial, hJperp_trial, h_l1_trial] = ...
                phaseI_violation_values(c_rest_trial, cl, cu, idxJ, idxJperp);
            [h_rest_trial, ~] = constraint_violation(c_rest_trial, cl, cu);

            best_iter_x = x_rest_trial;
            best_iter_h = h_rest_trial;
            best_iter_h_l1 = h_l1_trial;
            pred_l1_trial = sum(max(ccanon_rest + Acanon_rest * drest, 0));
            true_vs_pred = h_l1_trial / max(pred_l1_trial, 1e-12);
            true_decrease = h_l1_trial < h_l1_now;
            pred_decrease = pred_l1_trial < h_l1_now;

            filter_ok = filter_accept(rest_filter_hJ, rest_filter_hJperp, hJ_trial, hJperp_trial);
            l1_decrease = h_l1_trial < h_l1_now;
            max_not_worse = h_rest_trial <= h_rest_now;
            accepted_rest = filter_ok && l1_decrease && max_not_worse;

            % Record the restoration trial point.
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
            % End restoration trial-point record.

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
                if h_rest_trial < 1e-8
                    stop_restoration = true;
                end
            else
                soc_accepted = false;
                max_soc_iter = 5;
                d_soc_prev = drest;
                h_soc_prev_l1 = h_l1_trial;
                for soc_iter = 1:max_soc_iter
                    [d_soc, socflag] = solve_restoration_soc_qp( ...
                        prob, x_rest, d_soc_prev, restfun, Acanon_rest, idxJperp, ...
                        cl, cu, lb_d_rest, ub_d_rest, options);

                    % Record the restoration SOC exit flag.
                    rest_soc_log_id = numel(log.restoration_soc) + 1;
                    log.restoration_soc(rest_soc_log_id).outer_iter = k;
                    log.restoration_soc(rest_soc_log_id).rest_iter = rest_iter;
                    log.restoration_soc(rest_soc_log_id).soc_iter = soc_iter;
                    log.restoration_soc(rest_soc_log_id).socflag = socflag;
                    log.restoration_soc(rest_soc_log_id).h_l1_before = h_soc_prev_l1;
                    log.restoration_soc(rest_soc_log_id).h_l1_after = NaN;
                    log.restoration_soc(rest_soc_log_id).soc_rate = NaN;
                    log.restoration_soc(rest_soc_log_id).accepted_soc = false;
                    % End restoration SOC exit-flag record.

                    if socflag == -2
                        break;
                    elseif socflag <= 0
                        log.failure_reason = "restoration_soc_qp_failed";
                        break;
                    end
                    x_soc_trial = x_rest + d_soc;
                    c_soc_trial = eval_cons(prob, 1:m, x_soc_trial);
                    [hJ_soc, hJperp_soc, h_l1_soc] = ...
                        phaseI_violation_values(c_soc_trial, cl, cu, idxJ, idxJperp);
                    [h_rest_soc, ~] = constraint_violation(c_soc_trial, cl, cu);

                    % Record the restoration SOC trial point.
                    log.restoration_soc(rest_soc_log_id).hJ_soc = hJ_soc;
                    log.restoration_soc(rest_soc_log_id).hJperp_soc = hJperp_soc;
                    log.restoration_soc(rest_soc_log_id).h_l1_after = h_l1_soc;
                    if h_soc_prev_l1 > restoration_tol
                        log.restoration_soc(rest_soc_log_id).soc_rate = h_l1_soc / h_soc_prev_l1;
                    end
                    if exist('d_soc', 'var') && ~isempty(d_soc)
                        log.restoration_soc(rest_soc_log_id).norm_d_soc_inf = norm(d_soc, inf);
                    end
                    % End restoration SOC trial-point record.

                    if h_rest_soc < best_iter_h
                        best_iter_h = h_rest_soc;
                        best_iter_h_l1 = h_l1_soc;
                        best_iter_x = x_soc_trial;
                    end
                    accepted_soc = filter_accept(rest_filter_hJ, rest_filter_hJperp, hJ_soc, hJperp_soc);
                    log.restoration_soc(rest_soc_log_id).accepted_soc = accepted_soc;
                    if accepted_soc
                        x_rest = x_soc_trial;
                        restoration_success = true;
                        soc_accepted = true;
                        [rest_filter_hJ, rest_filter_hJperp] = filter_update( ...
                            rest_filter_hJ, rest_filter_hJperp, hJ_soc, hJperp_soc);
                        hJ_trial = hJ_soc;
                        hJperp_trial = hJperp_soc;
                        h_l1_trial = h_l1_soc;
                        if h_l1_trial < rest_upper_bound
                            rest_upper_bound = max(h_l1_trial, rest_upper_bound / 10);
                        end
                        if ~last_main_rejected && norm(d_soc,inf) >= 0.8 * rho
                            rho = min(2*rho, 50);
                        end
                        if h_l1_trial < 1e-8
                            stop_restoration = true;
                        end
                        break;
                    end
                    if h_l1_soc < 1e-8
                        break;
                    end
                    if h_soc_prev_l1 > restoration_tol
                        soc_rate = h_l1_soc / h_soc_prev_l1;
                        if soc_rate > 0.25
                            break;
                        end
                    end
                    d_soc_prev = d_soc;
                    h_soc_prev_l1 = h_l1_soc;
                end
                if ~soc_accepted
                    if J_changed && best_iter_h_l1 < rest_upper_bound
                        x_rest = best_iter_x;
                        restoration_success = true;
                        rest_filter_hJ = [];
                        rest_filter_hJperp = [];
                        rest_upper_bound = max(best_iter_h_l1, rest_upper_bound / 10);
                        if best_iter_h_l1 < 1e-8
                            stop_restoration = true;
                        end
                    else
                        rho = max(0.5 * rho, rho_min);
                    end
                end
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
                warning('Restoration loop failed at iteration %d', k);
                break;
            end
        end
    end
    x_trial = x_current + d;
    f_trial = eval_obj(prob, x_trial);
    c_trial = eval_cons(prob, 1:m, x_trial);
    [h_current, viol_current] = constraint_violation(cval, cl, cu);
    [h_trial, viol_trial] = constraint_violation(c_trial, cl, cu);

    % Record the main trial point.
    log.iter(k).f_trial = f_trial;
    log.iter(k).h_trial = h_trial;
    if exist('d', 'var') && ~isempty(d)
        log.iter(k).norm_d_inf = norm(d, inf);
    end
    % End main trial-point record.

    if restoration_phase
        accepted = h_trial < h_current;
    else
        accepted = h_trial <= filter_upper_bound && filter_accept(filter_f, filter_h, f_trial, h_trial);

        if first_after_restoration && h_trial > h_current
            accepted = false;
        end
    end
    best_soc_available = false;
    best_soc_d = [];
    best_soc_lambda = lambda;
    best_soc_x = [];
    best_soc_f = [];
    best_soc_c = [];
    best_soc_h = Inf;
    best_soc_merit = Inf;
    if ~accepted && ~restoration_phase && h_trial > 0
        max_main_soc_iter = 20;
        d_soc_prev = d;
        h_soc_prev = h_trial;
        soc_iter = 0;
        mu_soc = max(norm(lambda_before, inf), 1e-6);
        mu_soc = 10^ceil(log10(mu_soc));
        mu_soc = min(max(mu_soc, 1e-6), 1e6);
        while true
            soc_iter = soc_iter + 1;
            c_soc_prev = eval_cons(prob, 1:m, x_current + d_soc_prev);
            ccanon_soc_prev = canonical_value(c_soc_prev, cl, cu);
            Asoc = Acanon;
            bsoc = -ccanon_soc_prev + Acanon * d_soc_prev;
            [d_soc, ~, socflag, ~, lambda_soc] = fmincon( ...
                qpfun, d_soc_prev, Asoc, bsoc, [], [], ...
                lb_d, ub_d, [], options);

            % Record the main SOC exit flag.
            main_soc_log_id = numel(log.main_soc) + 1;
            log.main_soc(main_soc_log_id).outer_iter = k;
            log.main_soc(main_soc_log_id).soc_iter = soc_iter;
            log.main_soc(main_soc_log_id).socflag = socflag;
            log.main_soc(main_soc_log_id).h_before = h_soc_prev;
            log.main_soc(main_soc_log_id).h_after = NaN;
            log.main_soc(main_soc_log_id).f_soc_trial = NaN;
            log.main_soc(main_soc_log_id).soc_merit = NaN;
            log.main_soc(main_soc_log_id).norm_d_soc_inf = NaN;
            log.main_soc(main_soc_log_id).soc_rate = NaN;
            log.main_soc(main_soc_log_id).accepted_soc = false;
            % End main SOC exit-flag record.

            if socflag <= 0
                break;
            end
            x_soc_trial = x_current + d_soc;
            f_soc_trial = eval_obj(prob, x_soc_trial);
            c_soc_trial = eval_cons(prob, 1:m, x_soc_trial);
            [h_soc_trial, ~] = constraint_violation(c_soc_trial, cl, cu);
            soc_merit = f_soc_trial + mu_soc * h_soc_trial;

            % Record the main SOC trial point.
            log.main_soc(main_soc_log_id).h_after = h_soc_trial;
            log.main_soc(main_soc_log_id).f_soc_trial = f_soc_trial;
            log.main_soc(main_soc_log_id).soc_merit = soc_merit;
            log.main_soc(main_soc_log_id).norm_d_soc_inf = norm(d_soc, inf);
            if h_soc_prev > 1e-8
                log.main_soc(main_soc_log_id).soc_rate = h_soc_trial / h_soc_prev;
            else
                log.main_soc(main_soc_log_id).soc_rate = NaN;
            end
            % End main SOC trial-point record.

            if h_soc_trial <= h_soc_prev && soc_merit < best_soc_merit
                best_soc_available = true;
                best_soc_d = d_soc;
                best_soc_lambda = lambda_soc;
                best_soc_x = x_soc_trial;
                best_soc_f = f_soc_trial;
                best_soc_c = c_soc_trial;
                best_soc_h = h_soc_trial;
                best_soc_merit = soc_merit;
            end
          
            accepted_soc = h_soc_trial <= filter_upper_bound ...
                && h_soc_trial <= h_soc_prev ...
                && filter_accept(filter_f, filter_h, f_soc_trial, h_soc_trial);
            log.main_soc(main_soc_log_id).accepted_soc = accepted_soc;
            if accepted_soc
                d = d_soc;
                lambda = lambda_soc;
                x_trial = x_soc_trial;
                f_trial = f_soc_trial;
                c_trial = c_soc_trial;
                h_trial = h_soc_trial;
                accepted = true;
                break;
            end
            if h_soc_trial < 1e-8
                break;
            end
            if soc_iter >= max_main_soc_iter
                break;
            end
            if h_soc_prev > 1e-8
                soc_rate = h_soc_trial / h_soc_prev;
                if soc_rate > 0.25
                    break;
                end
            end
            d_soc_prev = d_soc;
            h_soc_prev = h_soc_trial;
        end
    end
    if ~accepted && ~restoration_phase && first_after_restoration ...
            && best_soc_available && best_soc_h <= h_current
        d = best_soc_d;
        lambda = best_soc_lambda;
        x_trial = best_soc_x;
        f_trial = best_soc_f;
        c_trial = best_soc_c;
        h_trial = best_soc_h;
        accepted = true;
        filter_f = [];
        filter_h = [];
        filter_upper_bound = max(best_soc_h, filter_upper_bound / 10);
    end
    if accepted
        log.iter(k).accepted = true;
        x = x_trial;

        last_main_rejected = false;
        if ~restoration_phase && isfield(lambda,'ineqlin') && length(lambda.ineqlin) >= m
            lambda_current = canonical_lambda_to_original(lambda.ineqlin, cl, cu, m);
        end
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
% End final-result record.
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
result.num_restoration_soc = numel(log.restoration_soc);
result.num_main_soc = numel(log.main_soc);
result.status = classify_run_status(log.final.success, h_final, log.failure_reason, cfg);
result.failure_reason = log.failure_reason;
result.failed_rest_message = "";
result.failed_restflag = NaN;
result.failed_rest_constrviolation = NaN;
result.failed_rest_firstorderopt = NaN;

if ~isempty(log.restoration)
    failed_rest_idx = find([log.restoration.restflag] <= 0, 1, 'last');
    if ~isempty(failed_rest_idx)
        result.failed_rest_message = log.restoration(failed_rest_idx).rest_message;
        result.failed_restflag = log.restoration(failed_rest_idx).restflag;
        result.failed_rest_constrviolation = log.restoration(failed_rest_idx).rest_constrviolation;
        result.failed_rest_firstorderopt = log.restoration(failed_rest_idx).rest_firstorderopt;
    end
end

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

function lambda_original = canonical_lambda_to_original(lambda_ineq, cl, cu, m)
lambda_original = zeros(m, 1);
idxU = find(isfinite(cu));
idxL = find(isfinite(cl));
nU = numel(idxU);
nL = numel(idxL);

if numel(lambda_ineq) >= nU
    lambda_original(idxU) = lambda_original(idxU) + lambda_ineq(1:nU);
end

if numel(lambda_ineq) >= nU + nL
    lambda_original(idxL) = lambda_original(idxL) - lambda_ineq(nU+1:nU+nL);
end
end

function [f, grad] = quad_objective(d, H, g)
f = 0.5 * d' * H * d + g' * d;
if nargout > 1
    grad = H * d + g;
end
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
    solve_elastic_lp_partition(Acanon, ccanon, lb_d, ub_d, n, restoration_tol, ...
    lp_options)

ncanon = length(ccanon);

felastic = [zeros(n,1); ones(ncanon,1)];

Aelastic = [Acanon, -eye(ncanon)];
belastic = -ccanon;

lbelastic = [lb_d; zeros(ncanon,1)];
ubelastic = [ub_d; Inf(ncanon,1)];

elasticfval = NaN;
elasticoutput = struct();

[zelastic, elasticfval, elasticflag, elasticoutput] = linprog( ...
    felastic, Aelastic, belastic, [], [], ...
    lbelastic, ubelastic, lp_options);

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


function [drest, restflag, restfun, restfval, restoutput] = solve_phase1_qp( ...
    prob, x_rest, lambda_before, Acanon_rest, ccanon_rest, ...
    idxJ, idxJperp, lb_d_rest, ub_d_rest, n, ...
    rest_hessian_mode, options)

switch rest_hessian_mode
    case "identity"
        Hrest = eye(n);

    case "lagrangian"
        Hrest = eval_hesslag(prob, lambda_before, x_rest);
        Hrest = 0.5 * (Hrest + Hrest');

    case "positive_lagrangian"
        Hrest = eval_hesslag(prob, lambda_before, x_rest);
        Hrest = 0.5 * (Hrest + Hrest');

        eigHrest_raw = eig(Hrest);
        minEigHrest_raw = min(eigHrest_raw);
        if minEigHrest_raw < 1e-6
            Hrest = Hrest + (-minEigHrest_raw + 1e-6) * eye(n);
        end

    otherwise
        error('Unknown rest_hessian_mode: %s', rest_hessian_mode);
end

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

restfun = @(drest) quad_objective(drest, Hrest, frest);
drest0 = zeros(n,1);

[drest, restfval, restflag, restoutput] = fmincon( ...
    restfun, drest0, Arest, brest, [], [], ...
    lb_d_rest, ub_d_rest, [], options);

end

function [hJ, hJperp, h_l1] = phaseI_violation_values(c, cl, cu, idxJ, idxJperp)
ccanon = canonical_value(c, cl, cu);

hJ = sum(max(ccanon(idxJ), 0));
hJperp = sum(max(ccanon(idxJperp), 0));
h_l1 = sum(max(ccanon, 0));
end

function [d_soc, socflag] = solve_restoration_soc_qp( ...
    prob, x_rest, d_soc_prev, restfun, Acanon_rest, idxJperp, ...
    cl, cu, lb_d_rest, ub_d_rest, options)

c_soc_prev = eval_cons(prob, 1:length(cl), x_rest + d_soc_prev);
ccanon_soc_prev = canonical_value(c_soc_prev, cl, cu);

Arest_soc = [];
brest_soc = [];

if any(idxJperp)
    Arest_soc = Acanon_rest(idxJperp,:);
    brest_soc = -ccanon_soc_prev(idxJperp) ...
        + Acanon_rest(idxJperp,:) * d_soc_prev;
end

[d_soc, ~, socflag, ~] = fmincon( ...
    restfun, d_soc_prev, Arest_soc, brest_soc, [], [], ...
    lb_d_rest, ub_d_rest, [], options);

end
