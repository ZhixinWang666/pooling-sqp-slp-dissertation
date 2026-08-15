%% Batch lockstep trace for QN-SQP: fmincon-sqp vs quadprog


problem_name = "fould3";
nlfile = "foulds3.nl";
saved_starts_file = 'fould3_starting_point';
quad_algorithm = "active-set";  
max_trace_iter = 80;
step_diff_tol = 1e-3;
lambda_diff_tol = 1;
B_diff_tol = Inf;
verbose = false;
save_trace = false;
print_all_runs = false;

S = load(saved_starts_file, 'starts');
num_runs = size(S.starts, 2);
run_ids = 1:num_runs;

branch_summary = repmat(make_empty_branch_summary(), num_runs, 1);

fprintf('\n=== Batch lockstep trace | %s | quadprog=%s | runs=%d ===\n', ...
    problem_name, quad_algorithm, num_runs);
if print_all_runs
    fprintf('%5s %8s %35s %12s %12s %12s %8s %8s %8s %8s\n', ...
        'run', 'iter', 'reason', 'last_d', 'last_lam', 'last_B', ...
        'acc_f', 'acc_q', 'soc_f', 'soc_q');
end

for ii = 1:num_runs
    run_id = run_ids(ii);

    run('Diagnose_QN_SQP_LockstepTrace.m');

    last_row = trace(end);
    branch_summary(ii).run_id = run_id;
    branch_summary(ii).divergence_iter = divergence_iter;
    branch_summary(ii).divergence_reason = divergence_reason;
    branch_summary(ii).last_norm_d_diff_inf = last_row.norm_d_diff_inf;
    branch_summary(ii).last_norm_lambda_diff_inf = last_row.norm_lambda_diff_inf;
    branch_summary(ii).last_norm_B_diff_inf = last_row.norm_B_diff_inf;
    branch_summary(ii).fmincon_accepted = last_row.fmincon_accepted;
    branch_summary(ii).quadprog_accepted = last_row.quadprog_accepted;
    branch_summary(ii).fmincon_would_enter_main_soc = last_row.fmincon_would_enter_main_soc;
    branch_summary(ii).quadprog_would_enter_main_soc = last_row.quadprog_would_enter_main_soc;

    if print_all_runs
        fprintf('%5d %8g %35s %12.3e %12.3e %12.3e %8d %8d %8d %8d\n', ...
            branch_summary(ii).run_id, branch_summary(ii).divergence_iter, ...
            branch_summary(ii).divergence_reason, ...
            branch_summary(ii).last_norm_d_diff_inf, ...
            branch_summary(ii).last_norm_lambda_diff_inf, ...
            branch_summary(ii).last_norm_B_diff_inf, ...
            branch_summary(ii).fmincon_accepted, ...
            branch_summary(ii).quadprog_accepted, ...
            branch_summary(ii).fmincon_would_enter_main_soc, ...
            branch_summary(ii).quadprog_would_enter_main_soc);
    end
end

reasons = string({branch_summary.divergence_reason});
unique_reasons = unique(reasons);
fprintf('\n=== Branch reason counts ===\n');
for i = 1:numel(unique_reasons)
    fprintf('%35s : %d\n', unique_reasons(i), sum(reasons == unique_reasons(i)));
end

fprintf('\n=== Representative branch examples ===\n');
for i = 1:numel(unique_reasons)
    idx = find(reasons == unique_reasons(i), 1, 'first');
    row = branch_summary(idx);
    fprintf('run %d, iter %g, reason=%s, d_diff=%.3e, lambda_diff=%.3e, B_diff=%.3e, acc_f/acc_q=%d/%d, soc_f/soc_q=%d/%d\n', ...
        row.run_id, row.divergence_iter, row.divergence_reason, ...
        row.last_norm_d_diff_inf, row.last_norm_lambda_diff_inf, ...
        row.last_norm_B_diff_inf, row.fmincon_accepted, row.quadprog_accepted, ...
        row.fmincon_would_enter_main_soc, row.quadprog_would_enter_main_soc);
end

fprintf('\nInterpretation: first meaningful differences occur later in the QN-SQP path, through step/lambda differences, filter/SOC branching, or restoration-related branching.\n');

timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
result_file = char("QN_SQP_lockstep_trace_all_" + problem_name + "_" + quad_algorithm + "_" + timestamp + ".mat");
save(result_file, 'branch_summary', 'problem_name', 'nlfile', 'saved_starts_file', ...
    'quad_algorithm', 'max_trace_iter', 'step_diff_tol', 'lambda_diff_tol', 'B_diff_tol');
fprintf('Saved batch lockstep summary to %s\n', result_file);


function row = make_empty_branch_summary()
row = struct();
row.run_id = NaN;
row.divergence_iter = NaN;
row.divergence_reason = "";
row.last_norm_d_diff_inf = NaN;
row.last_norm_lambda_diff_inf = NaN;
row.last_norm_B_diff_inf = NaN;
row.fmincon_accepted = false;
row.quadprog_accepted = false;
row.fmincon_would_enter_main_soc = false;
row.quadprog_would_enter_main_soc = false;
end
