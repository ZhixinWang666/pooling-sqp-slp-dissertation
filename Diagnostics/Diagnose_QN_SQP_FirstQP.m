%% Diagnose the first QP in QN-SQP: fmincon vs quadprog

rng('default');

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
if size(starts, 1) ~= n
    error('Saved starts dimension does not match this problem.');
end
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

diagnostics = repmat(make_empty_diag(), numel(run_ids), 1);

fprintf('\n=== First QP diagnostic: QN-SQP | %s | %s ===\n', problem_name, nlfile);
if print_all_rows
    fprintf('%5s %7s %7s %12s %12s %12s %12s %12s %12s %12s\n', ...
        'run', 'flag_f', 'flag_q', 'normdiff', 'qobj_f', 'qobj_q', ...
        'qobj_gap', 'viol_f', 'viol_q', 'mineigH');
end

for ii = 1:numel(run_ids)
    run_id = run_ids(ii);
    x = starts(:, run_id);
    diag = compare_first_qp(prob, n, m, bl, bu, cl, cu, x, ...
        fmincon_options, quadprog_options);
    diag.run_id = run_id;
    diagnostics(ii) = diag;

    if print_all_rows
        fprintf('%5d %7d %7d %12.3e %12.6g %12.6g %12.3e %12.3e %12.3e %12.3e\n', ...
            diag.run_id, diag.fmincon_exitflag, diag.quadprog_exitflag, ...
            diag.norm_d_diff_inf, diag.qobj_fmincon, diag.qobj_quadprog, ...
            diag.qobj_gap_fminusq, diag.maxviol_fmincon, diag.maxviol_quadprog, ...
            diag.min_eig_H);
    end
end

fprintf('\n=== First QP max diagnostics over %d starts ===\n', numel(diagnostics));
all_flags_ok = all([diagnostics.fmincon_exitflag] == 1) && all([diagnostics.quadprog_exitflag] == 1);
[max_step_diff, worst_step_idx] = max([diagnostics.norm_d_diff_inf]);
[max_qobj_gap, worst_qobj_idx] = max(abs([diagnostics.qobj_gap_fminusq]));
max_viol_f = max([diagnostics.maxviol_fmincon]);
max_viol_q = max([diagnostics.maxviol_quadprog]);
fprintf('All first QPs solved successfully by both solvers: %d\n', all_flags_ok);
fprintf('max ||d_fmincon - d_quadprog||_inf = %.3e (run %d)\n', ...
    max_step_diff, diagnostics(worst_step_idx).run_id);
fprintf('max abs(qobj_fmincon - qobj_quadprog) = %.3e (run %d)\n', ...
    max_qobj_gap, diagnostics(worst_qobj_idx).run_id);
fprintf('max linear violation = %.3e (fmincon), %.3e (quadprog)\n', ...
    max_viol_f, max_viol_q);
fprintf('min eigenvalue range of H = [%.3e, %.3e]\n', ...
    min([diagnostics.min_eig_H]), max([diagnostics.min_eig_H]));
fprintf('Interpretation: the first convex QP solutions agree within numerical tolerance.\n');

timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
result_file = char("QN_SQP_firstQP_diagnostic_" + problem_name + "_" + timestamp + ".mat");
save(result_file, 'diagnostics', 'problem_name', 'nlfile', 'saved_starts_file', 'run_ids');
fprintf('Saved diagnostics to %s\n', result_file);


function diag = compare_first_qp(prob, n, m, bl, bu, cl, cu, x, fmincon_options, quadprog_options)
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

[d_fmincon, qobj_fmincon, fmincon_exitflag, fmincon_output, fmincon_lambda] = fmincon( ...
    qpfun, d0, Aqp, bqp, [], [], lb_d, ub_d, [], fmincon_options);

[d_quadprog, qobj_quadprog, quadprog_exitflag, quadprog_output, quadprog_lambda] = quadprog( ...
    H, g(:), Aqp, bqp, [], [], lb_d, ub_d, d0, quadprog_options);

diag = make_empty_diag();
diag.min_eig_H = min(eig(H));
diag.norm_d_diff_inf = norm(d_fmincon - d_quadprog, inf);
diag.norm_d_diff_2 = norm(d_fmincon - d_quadprog, 2);
diag.qobj_fmincon = qobj_fmincon;
diag.qobj_quadprog = qobj_quadprog;
diag.qobj_gap_fminusq = qobj_fmincon - qobj_quadprog;
diag.maxviol_fmincon = max_linear_qp_violation(Aqp, bqp, lb_d, ub_d, d_fmincon);
diag.maxviol_quadprog = max_linear_qp_violation(Aqp, bqp, lb_d, ub_d, d_quadprog);
diag.pred_canon_max_fmincon = max(ccanon + Acanon * d_fmincon);
diag.pred_canon_max_quadprog = max(ccanon + Acanon * d_quadprog);
diag.fmincon_exitflag = fmincon_exitflag;
diag.quadprog_exitflag = quadprog_exitflag;
diag.fmincon_iterations = get_output_field(fmincon_output, 'iterations');
diag.quadprog_iterations = get_output_field(quadprog_output, 'iterations');
diag.fmincon_firstorderopt = get_output_field(fmincon_output, 'firstorderopt');
diag.quadprog_firstorderopt = get_output_field(quadprog_output, 'firstorderopt');
diag.fmincon_message = get_output_message(fmincon_output);
diag.quadprog_message = get_output_message(quadprog_output);
diag.d_fmincon = d_fmincon;
diag.d_quadprog = d_quadprog;
diag.fmincon_lambda = fmincon_lambda;
diag.quadprog_lambda = quadprog_lambda;
end


function diag = make_empty_diag()
diag = struct();
diag.run_id = NaN;
diag.min_eig_H = NaN;
diag.norm_d_diff_inf = NaN;
diag.norm_d_diff_2 = NaN;
diag.qobj_fmincon = NaN;
diag.qobj_quadprog = NaN;
diag.qobj_gap_fminusq = NaN;
diag.maxviol_fmincon = NaN;
diag.maxviol_quadprog = NaN;
diag.pred_canon_max_fmincon = NaN;
diag.pred_canon_max_quadprog = NaN;
diag.fmincon_exitflag = NaN;
diag.quadprog_exitflag = NaN;
diag.fmincon_iterations = NaN;
diag.quadprog_iterations = NaN;
diag.fmincon_firstorderopt = NaN;
diag.quadprog_firstorderopt = NaN;
diag.fmincon_message = "";
diag.quadprog_message = "";
diag.d_fmincon = [];
diag.d_quadprog = [];
diag.fmincon_lambda = struct();
diag.quadprog_lambda = struct();
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


function maxviol = max_linear_qp_violation(A, b, lb, ub, d)
viol = [];
if ~isempty(A)
    viol = [viol; A * d - b];
end
viol = [viol; lb - d; d - ub];
viol(~isfinite(viol)) = 0;
maxviol = max(max(viol, 0));
end


function value = get_output_field(output, name)
if isfield(output, name) && ~isempty(output.(name))
    value = output.(name);
else
    value = NaN;
end
end


function message = get_output_message(output)
if isfield(output, 'message') && ~isempty(output.message)
    message = string(output.message);
else
    message = "";
end
end
