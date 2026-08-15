
rng();

target_f = (sqrt(5) - 1)^2;
cfg = make_problem("Simpleproblem", "Simpleproblem.nl", [], "manual", target_f);
cfg.x_star = [2; 1] / sqrt(5);

prob = NonlinearProblem(char(cfg.nlfile));
n = get_nvar(prob);
m = get_ncon(prob);
bl = get_bl(prob);
bu = get_bu(prob);
cl = get_cl(prob);
cu = get_cu(prob);

x = [0.6; 0.2];

run_timer = tic;
[result, log] = solve_pure_slp_once(prob, n, m, bl, bu, cl, cu, x, cfg);
result.solve_time = toc(run_timer);

fprintf('\n=== Pure SLP validation | %s | %s | mode=%s ===\n', ...
    char(cfg.name), char(cfg.nlfile), char(cfg.start_mode));
fprintf('success       : %d\n', result.success);
fprintf('status        : %s\n', char(result.status));
fprintf('reason        : %s\n', char(result.failure_reason));
fprintf('final f       : %.12g\n', result.final_f);
fprintf('final h       : %.3e\n', result.final_h);
fprintf('final rho     : %.3e\n', result.final_rho);
fprintf('iterations    : %d\n', result.iterations);
fprintf('solve time    : %.4f s\n', result.solve_time);

out_dir = "validation_figures";
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
file_prefix = lower(regexprep(char(cfg.name), '[^A-Za-z0-9]+', '_'));
save(fullfile(out_dir, sprintf('%s_pureslp_validation.mat', file_prefix)), ...
    'result', 'log', 'x', 'cfg');

plot_validation(log, result, cfg);


function plot_validation(log, result, cfg)
if isempty(log.iter)
    warning('No iteration history was recorded; validation plots were not created.');
    return
end

iter = log.iter;
k = [iter.k];
f = [iter.f];
h = [iter.h_current];
rho = [iter.rho];
norm_d = [iter.norm_d_inf];
accepted = [iter.accepted];

h_plot = max(h, 1e-16);
tr_activity = norm_d ./ rho;
method_label = 'Pure SLP';

out_dir = "validation_figures";
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
file_prefix = lower(regexprep(char(cfg.name), '[^A-Za-z0-9]+', '_'));

tailN = min(20, numel(k));
tail_idx = (numel(k) - tailN + 1):numel(k);

fig3 = figure('Name', 'Pure SLP trust-region activity');
tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
semilogy(k(tail_idx), max(rho(tail_idx), 1e-16), '-o', ...
    'LineWidth', 1.2, 'MarkerSize', 4);
hold on;
semilogy(k(tail_idx), max(norm_d(tail_idx), 1e-16), '-s', ...
    'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Value');
legend('\Delta_k', '||d_k||_\infty', 'Location', 'best');
grid on;

nexttile;
plot(k(tail_idx), tr_activity(tail_idx), '-^', ...
    'Color', [0.9290 0.6940 0.1250], ...
    'LineWidth', 1.2, 'MarkerSize', 4);
xlabel('Iteration k');
ylabel('||d_k||_\infty / \Delta_k');
ylim([0, 1.1 * max(1, max(tr_activity(tail_idx)))]);
grid on;
export_validation_figure(fig3, out_dir, sprintf('%s_pureslp_trust_region_activity', file_prefix));

if isfield(cfg, 'x_star')
    x_star = cfg.x_star(:);
    x_hist = reshape([iter.x_current], numel(x_star), []);

    err = vecnorm(x_hist - x_star, 2, 1);
    err_plot = max(err, 1e-16);
    k_err = k;

    fig4 = figure('Name', 'Pure SLP exact-solution error');
    semilogy(k_err, err_plot, '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
    xlabel('Iteration k');
    ylabel('||x_k - x_*||_2');
    grid on;
    export_validation_figure(fig4, out_dir, sprintf('%s_pureslp_error_history', file_prefix));

    if numel(err) >= 3
        err_k = err(1:end-1);
        err_next = err(2:end);
        linear_ratio = err_next ./ max(err_k, 1e-16);
        quadratic_ratio = err_next ./ max(err_k.^2, 1e-32);
        Tconv = table(k_err(1:end-1)', err_k(:), err_next(:), ...
            linear_ratio(:), quadratic_ratio(:), ...
            'VariableNames', {'k','error_k','error_next','linear_ratio','quadratic_ratio'});
        fprintf('\n=== exact solution error and convergence ratios ===\n');
        disp(Tconv(max(1, height(Tconv)-9):height(Tconv), :));
    end
end

fprintf('\n=== validation plots saved ===\n');
fprintf('%s\n', fullfile(out_dir, sprintf('%s_pureslp_trust_region_activity.png', file_prefix)));
fprintf('%s\n', fullfile(out_dir, sprintf('%s_pureslp_trust_region_activity.pdf', file_prefix)));
if isfield(cfg, 'x_star')
    fprintf('%s\n', fullfile(out_dir, sprintf('%s_pureslp_error_history.png', file_prefix)));
    fprintf('%s\n', fullfile(out_dir, sprintf('%s_pureslp_error_history.pdf', file_prefix)));
end

fprintf('\n=== selected trust-region diagnostics ===\n');
diag_table = table(k(:), f(:), h(:), rho(:), norm_d(:), tr_activity(:), accepted(:), ...
    'VariableNames', {'k','f','h','rho','norm_d','norm_d_over_rho','accepted'});
disp(diag_table(max(1, height(diag_table)-9):height(diag_table), :));
end


function export_validation_figure(fig, out_dir, basename)
set(fig, 'Color', 'w');
exportgraphics(fig, fullfile(out_dir, [char(basename), '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(out_dir, [char(basename), '.pdf']), 'ContentType', 'vector');
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


function [result, log] = solve_pure_slp_once(prob, n, m, bl, bu, cl, cu, x, cfg)
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
log.final = struct();
log.failure_reason = "";

rho = cfg.rho0;
rho_min = 1e-8;
max_iter = cfg.max_iter;

filter_f = [];
filter_h = [];
c_init = eval_cons(prob, 1:m, x);
[h_init, ~] = constraint_violation(c_init, cl, cu);
filter_upper_bound = max(0.5, 1.25 * h_init);
lp_options = optimoptions('linprog','Display','none');

small_progress_count = 0;
small_progress_limit = 10;

for k = 1:max_iter
    x_current = x;
    f = eval_obj(prob, x_current);
    g = eval_obj_grad(prob, x_current);
    cval = eval_cons(prob, 1:m, x_current);
    J = eval_jac_val(prob, 1:m, x_current);
    J = J';
    [Acanon, ccanon] = canonical_linearization(cval, J, cl, cu);
    [h_current, ~] = constraint_violation(cval, cl, cu);

    log.iter(k).k = k;
    log.iter(k).f = f;
    log.iter(k).h_current = h_current;
    log.iter(k).rho = rho;
    log.iter(k).norm_g_inf = norm(g, inf);
    log.iter(k).main_qp_exitflag = NaN;
    log.iter(k).norm_d_inf = NaN;
    log.iter(k).accepted = false;
    log.iter(k).f_trial = NaN;
    log.iter(k).h_trial = NaN;
    log.iter(k).pred_obj_change = NaN;
    log.iter(k).actual_obj_change = NaN;
    log.iter(k).x_current = x_current;
    log.iter(k).x_trial = NaN(n,1);

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

    [d, ~, exitflag] = linprog(g(:), Aqp, bqp, [], [], lb_d, ub_d, lp_options);
    log.iter(k).main_qp_exitflag = exitflag;

    if exitflag <= 0 || isempty(d)
        rho = max(0.5 * rho, rho_min);
        if rho <= rho_min
            log.failure_reason = "lp_failed";
            break;
        end
        continue;
    end

    log.iter(k).norm_d_inf = norm(d, inf);

    x_trial = x_current + d;
    f_trial = eval_obj(prob, x_trial);
    c_trial = eval_cons(prob, 1:m, x_trial);
    [h_trial, ~] = constraint_violation(c_trial, cl, cu);

    log.iter(k).f_trial = f_trial;
    log.iter(k).h_trial = h_trial;
    log.iter(k).pred_obj_change = g' * d;
    log.iter(k).actual_obj_change = f_trial - f;
    log.iter(k).x_trial = x_trial;

    accepted = h_trial <= filter_upper_bound && ...
        filter_accept(filter_f, filter_h, f_trial, h_trial);

    if accepted
        log.iter(k).accepted = true;
        x = x_trial;
        [filter_f, filter_h] = filter_update(filter_f, filter_h, f_trial, h_trial);
    else
        rho = max(0.5 * rho, rho_min);
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
log.final.rho = rho;
if ~isempty(log.iter) && isfield(log.iter(end), 'norm_d_inf')
    log.final.norm_d_inf = log.iter(end).norm_d_inf;
else
    log.final.norm_d_inf = NaN;
end
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

result = struct();
result.success = log.final.success;
result.final_f = f_final;
result.final_h = h_final;
result.final_rho = log.final.rho;
result.final_norm_d = log.final.norm_d_inf;
result.max_viol = max_viol_final;
result.max_viol_idx = max_viol_final_idx;
result.iterations = k;
result.x0 = x0;
result.x_final = x;
result.f0 = f0;
result.h0 = h0;
result.status = classify_run_status(log.final.success, h_final, log.failure_reason, cfg);
result.failure_reason = log.failure_reason;
end


function print_single_result(result, log, prob, m, cl, cu)
disp('=== final objective ===')
disp(result.final_f)

disp('=== final feasibility ===')
disp(struct( ...
    'x', result.x_final, ...
    'f', result.final_f, ...
    'h', result.final_h, ...
    'rho', result.final_rho, ...
    'norm_d', result.final_norm_d, ...
    'max_viol', result.max_viol, ...
    'max_viol_idx', result.max_viol_idx, ...
    'iterations', result.iterations, ...
    'success', result.success, ...
    'status', result.status))

disp('=== failure reason ===')
disp(log.failure_reason)

c_final = eval_cons(prob, 1:m, result.x_final);
[~, viol_final] = constraint_violation(c_final, cl, cu);
[viol_sorted, idx_sorted] = sort(viol_final, 'descend');
topN = min(10, numel(viol_sorted));
top_table = zeros(topN, 5);
for ii = 1:topN
    idx = idx_sorted(ii);
    top_table(ii,:) = [idx, viol_sorted(ii), c_final(idx), cl(idx), cu(idx)];
end

disp('=== top violated constraints: [idx, violation, c, cl, cu] ===')
disp(top_table)

if ~isempty(log.iter)
    disp('=== last main iterations ===')
    iter_from = max(1, length(log.iter) - 14);
    fprintf('%5s %14s %12s %10s %8s %10s %8s %14s %12s %12s %12s\n', ...
        'k', 'f', 'h', 'rho', 'exit', 'norm_d', 'accept', ...
        'f_trial', 'h_trial', 'pred_df', 'act_df');
    for ii = iter_from:length(log.iter)
        fprintf('%5d %14.6g %12.3e %10.3e %8g %10.3e %8d %14.6g %12.3e %12.3e %12.3e\n', ...
            log.iter(ii).k, ...
            log.iter(ii).f, ...
            log.iter(ii).h_current, ...
            log.iter(ii).rho, ...
            log.iter(ii).main_qp_exitflag, ...
            log.iter(ii).norm_d_inf, ...
            log.iter(ii).accepted, ...
            log.iter(ii).f_trial, ...
            log.iter(ii).h_trial, ...
            log.iter(ii).pred_obj_change, ...
            log.iter(ii).actual_obj_change);
    end
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
