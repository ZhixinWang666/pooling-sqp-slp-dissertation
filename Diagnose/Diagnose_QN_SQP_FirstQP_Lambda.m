%% Diagnose first-QP multipliers in QN-SQP: fmincon vs quadprog


problem_name = "fould3";
nlfile = "foulds3.nl";
saved_starts_file = 'fould3_starting_point';
run_ids = [];
print_all_rows = false;

prob = NonlinearProblem(char(nlfile));
n = get_nvar(prob);
m = get_ncon(prob);
bl = get_bl(prob);
bu = get_bu(prob);
cl = get_cl(prob);
cu = get_cu(prob);

S = load(saved_starts_file, 'starts');
starts = S.starts;
if isempty(run_ids)
    run_ids = 1:size(starts, 2);
else
    run_ids = run_ids(run_ids >= 1 & run_ids <= size(starts, 2));
end

fmincon_options = optimoptions('fmincon', ...
    'Display','none', ...
    'Algorithm','sqp', ...
    'SpecifyObjectiveGradient',true, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-6, ...
    'ConstraintTolerance',1e-6);

quadprog_options = optimoptions('quadprog', ...
    'Display','none', ...
    'Algorithm','active-set', ...
    'MaxIterations',500, ...
    'OptimalityTolerance',1e-6, ...
    'ConstraintTolerance',1e-6);

lambda_diagnostics = repmat(make_empty_lambda_diag(), numel(run_ids), 1);

fprintf('\n=== First QP lambda diagnostic: QN-SQP | %s | %s ===\n', problem_name, nlfile);
if print_all_rows
    fprintf('%5s %12s %12s %12s %12s %12s\n', ...
        'run', 'norm_d_inf', 'lam_ineq_inf', 'lam_orig_inf', 'gradL_inf', 'active_gap');
end

for ii = 1:numel(run_ids)
    run_id = run_ids(ii);
    x = starts(:, run_id);
    diag = compare_first_qp_lambda(prob, n, m, bl, bu, cl, cu, x, ...
        fmincon_options, quadprog_options);
    diag.run_id = run_id;
    lambda_diagnostics(ii) = diag;

    if print_all_rows
        fprintf('%5d %12.3e %12.3e %12.3e %12.3e %12.3e\n', ...
            diag.run_id, diag.norm_d_diff_inf, diag.norm_lambda_ineq_diff_inf, ...
            diag.norm_lambda_original_diff_inf, diag.norm_gradL_new_diff_inf, ...
            diag.active_set_gap_inf);
    end
end

fprintf('\n=== First QP lambda max diagnostics over %d starts ===\n', numel(lambda_diagnostics));
[max_step_diff, worst_step_idx] = max([lambda_diagnostics.norm_d_diff_inf]);
[max_lam_ineq, worst_lam_ineq_idx] = max([lambda_diagnostics.norm_lambda_ineq_diff_inf]);
[max_lam_orig, worst_lam_orig_idx] = max([lambda_diagnostics.norm_lambda_original_diff_inf]);
[max_gradL, worst_gradL_idx] = max([lambda_diagnostics.norm_gradL_new_diff_inf]);
fprintf('max ||d_fmincon - d_quadprog||_inf = %.3e (run %d)\n', ...
    max_step_diff, lambda_diagnostics(worst_step_idx).run_id);
fprintf('max ||lambda_ineq_fmincon - lambda_ineq_quadprog||_inf = %.3e (run %d)\n', ...
    max_lam_ineq, lambda_diagnostics(worst_lam_ineq_idx).run_id);
fprintf('max ||lambda_original_fmincon - lambda_original_quadprog||_inf = %.3e (run %d)\n', ...
    max_lam_orig, lambda_diagnostics(worst_lam_orig_idx).run_id);
fprintf('max ||gradL_new_fmincon - gradL_new_quadprog||_inf = %.3e (run %d)\n', ...
    max_gradL, lambda_diagnostics(worst_gradL_idx).run_id);
fprintf('Interpretation: first-QP multipliers and the BFGS gradient input agree within numerical tolerance.\n');

timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
result_file = char("QN_SQP_firstQP_lambda_diagnostic_" + problem_name + "_" + timestamp + ".mat");
save(result_file, 'lambda_diagnostics', 'problem_name', 'nlfile', 'saved_starts_file', 'run_ids');
fprintf('Saved lambda diagnostics to %s\n', result_file);


function diag = compare_first_qp_lambda(prob, n, m, bl, bu, cl, cu, x, fmincon_options, quadprog_options)
B = eye(n);
rho = 1;

g = eval_obj_grad(prob, x);
cval = eval_cons(prob, 1:m, x);
J = eval_jac_val(prob, 1:m, x);
J = J';
H = 0.5 * (B + B');
[Acanon, ccanon] = canonical_linearization(cval, J, cl, cu);

Aqp = Acanon;
bqp = -ccanon;
lb_d = max(bl - x, -rho * ones(n,1));
ub_d = min(bu - x,  rho * ones(n,1));
d0 = zeros(n,1);
qpfun = @(d) quad_objective(d, H, g);

[d_f, ~, flag_f, ~, lambda_f] = fmincon( ...
    qpfun, d0, Aqp, bqp, [], [], lb_d, ub_d, [], fmincon_options);

[d_q, ~, flag_q, ~, lambda_q] = quadprog( ...
    H, g(:), Aqp, bqp, [], [], lb_d, ub_d, d0, quadprog_options);

lambda_ineq_f = get_lambda_ineq(lambda_f);
lambda_ineq_q = get_lambda_ineq(lambda_q);
lambda_original_f = canonical_lambda_to_original(lambda_ineq_f, cl, cu, m);
lambda_original_q = canonical_lambda_to_original(lambda_ineq_q, cl, cu, m);

x_trial_f = x + d_f;
x_trial_q = x + d_q;
g_new_f = eval_obj_grad(prob, x_trial_f);
g_new_q = eval_obj_grad(prob, x_trial_q);
J_new_f = eval_jac_val(prob, 1:m, x_trial_f);
J_new_q = eval_jac_val(prob, 1:m, x_trial_q);
J_new_f = J_new_f';
J_new_q = J_new_q';

gradL_new_f = g_new_f + J_new_f' * lambda_original_f;
gradL_new_q = g_new_q + J_new_q' * lambda_original_q;

active_res_f = Aqp * d_f - bqp;
active_res_q = Aqp * d_q - bqp;

diag = make_empty_lambda_diag();
diag.fmincon_exitflag = flag_f;
diag.quadprog_exitflag = flag_q;
diag.norm_d_diff_inf = norm(d_f - d_q, inf);
diag.norm_lambda_ineq_diff_inf = norm_with_padding(lambda_ineq_f, lambda_ineq_q, inf);
diag.norm_lambda_original_diff_inf = norm(lambda_original_f - lambda_original_q, inf);
diag.norm_gradL_new_diff_inf = norm(gradL_new_f - gradL_new_q, inf);
diag.active_set_gap_inf = norm(active_res_f - active_res_q, inf);
diag.lambda_ineq_fmincon = lambda_ineq_f;
diag.lambda_ineq_quadprog = lambda_ineq_q;
diag.lambda_original_fmincon = lambda_original_f;
diag.lambda_original_quadprog = lambda_original_q;
end


function diag = make_empty_lambda_diag()
diag = struct();
diag.run_id = NaN;
diag.fmincon_exitflag = NaN;
diag.quadprog_exitflag = NaN;
diag.norm_d_diff_inf = NaN;
diag.norm_lambda_ineq_diff_inf = NaN;
diag.norm_lambda_original_diff_inf = NaN;
diag.norm_gradL_new_diff_inf = NaN;
diag.active_set_gap_inf = NaN;
diag.lambda_ineq_fmincon = [];
diag.lambda_ineq_quadprog = [];
diag.lambda_original_fmincon = [];
diag.lambda_original_quadprog = [];
end


function lambda_ineq = get_lambda_ineq(lambda)
if isfield(lambda, 'ineqlin') && ~isempty(lambda.ineqlin)
    lambda_ineq = lambda.ineqlin(:);
else
    lambda_ineq = [];
end
end


function value = norm_with_padding(a, b, p)
na = numel(a);
nb = numel(b);
nmax = max(na, nb);
aa = zeros(nmax, 1);
bb = zeros(nmax, 1);
aa(1:na) = a(:);
bb(1:nb) = b(:);
value = norm(aa - bb, p);
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
f = 0.5 * d' * H * d + g(:)' * d;
if nargout > 1
    grad = H * d + g(:);
end
end


function [Acanon, ccanon] = canonical_linearization(c, J, cl, cu)
idxU = isfinite(cu);
idxL = isfinite(cl);
Acanon = [J(idxU,:); -J(idxL,:)];
ccanon = [c(idxU) - cu(idxU); cl(idxL) - c(idxL)];
end
