function x0 = generate_fould3_start(mode)


if nargin < 1
    mode = "low";
else
    mode = string(mode);
end

switch mode
    case "low"
        demand_scale = 0.10 + 0.20 * rand();
    case "medium"
        demand_scale = 0.35 + 0.45 * rand();
    case "boundary"
        demand_scale = 0.75 + 0.25 * rand();
    case "quality_hard"
        demand_scale = 0.45 + 0.45 * rand();
    otherwise
        demand_scale = 0.35 + 0.45 * rand();
end

num_c = 32;
num_o = 8;
num_p = 16;

conc_c = [
    1.0; 1.1; 1.2; 1.3
    1.1; 1.2; 1.3; 1.4
    1.2; 1.3; 1.4; 1.5
    1.3; 1.4; 1.5; 1.6
    1.4; 1.5; 1.6; 1.7
    1.5; 1.6; 1.7; 1.8
    1.6; 1.7; 1.8; 1.9
    1.7; 1.8; 1.9; 2.0
    ];
target_p = (1.05:0.05:1.80)';


pool_comps = {
    1:4
    5:8
    9:12
    13:16
    17:20
    21:24
    25:28
    29:32
    };


q_pool = [1.02; 1.16; 1.28; 1.40; 1.52; 1.64; 1.72; 1.78];

if mode == "quality_hard"
    q_pool = q_pool + 0.02 * randn(num_o,1);
else
    q_pool = q_pool + 0.01 * randn(num_o,1);
end

for o = 1:num_o
    comps = pool_comps{o};
    qmin = min(conc_c(comps));
    qmax = max(conc_c(comps));
    q_pool(o) = min(max(q_pool(o), qmin + 1e-3), qmax - 1e-3);
end


o2p = zeros(num_o, num_p);
for p = 1:num_p
    demand_p = demand_scale * (0.40 + 0.60 * rand());  % ub_p is 1
    feasible_pools = find(q_pool <= target_p(p) + 1e-10);

    if isempty(feasible_pools)
        [~, pool_id] = min(q_pool);
    else
        [~, order] = sort(q_pool(feasible_pools), 'descend');
        sorted_pools = feasible_pools(order);
        if numel(sorted_pools) >= 2 && rand() < 0.20
            pool_id = sorted_pools(2);
        else
            pool_id = sorted_pools(1);
        end
    end

    o2p(pool_id, p) = demand_p;
end

fmat = zeros(num_c, num_o);
C1_o = zeros(num_o,1);
for o = 1:num_o
    total_o = sum(o2p(o,:));
    comps = pool_comps{o};
    q_target = q_pool(o);

    f_pool = split_to_quality(total_o, q_target, conc_c(comps));
    fmat(comps, o) = f_pool;

    if total_o > 1e-12
        C1_o(o) = dot(conc_c(comps), f_pool) / sum(f_pool);
    else
        C1_o(o) = q_target;
    end
end

x_raw = sum(fmat, 2);

max_raw = max(x_raw);
if max_raw > 15
    factor = 15 / max_raw;
    fmat = factor * fmat;
    o2p = factor * o2p;
    x_raw = factor * x_raw;
    for o = 1:num_o
        total_o = sum(fmat(:,o));
        if total_o > 1e-12
            C1_o(o) = dot(conc_c, fmat(:,o)) / total_o;
        end
    end
end

f_col = [
    fmat(1:4,1)
    fmat(5:8,2)
    fmat(9:12,3)
    fmat(13:16,4)
    fmat(17:20,5)
    fmat(21:24,6)
    fmat(25:28,7)
    fmat(29:32,8)
    ];

o2p_col = reshape(o2p', [], 1);

x0 = [
    f_col
    o2p_col
    C1_o
    x_raw
    ];

assert(numel(x0) == 200, 'foulds3 start vector must have 200 variables.');


end

function f = split_to_quality(total_flow, q_target, comp_qualities)


num_comp = numel(comp_qualities);
f = zeros(num_comp,1);

if total_flow <= 1e-12
    return
end

[q_sorted, order] = sort(comp_qualities(:));
q_target = min(max(q_target, q_sorted(1)), q_sorted(end));

lower_pos = find(q_sorted <= q_target + 1e-12, 1, 'last');
upper_pos = find(q_sorted >= q_target - 1e-12, 1, 'first');

if lower_pos == upper_pos
    f(order(lower_pos)) = total_flow;
    return
end

q_low = q_sorted(lower_pos);
q_high = q_sorted(upper_pos);

if abs(q_high - q_low) <= 1e-12
    f(order(lower_pos)) = total_flow;
    return
end

weight_high = (q_target - q_low) / (q_high - q_low);
weight_low = 1 - weight_high;

f(order(lower_pos)) = total_flow * weight_low;
f(order(upper_pos)) = total_flow * weight_high;
end
