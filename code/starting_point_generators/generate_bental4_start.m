function x0 = generate_bental4_start(mode)


if nargin < 1
    mode = "low";
else
    mode = string(mode);
end

switch mode
    case "low"
        scale = 0.05 + 0.15 * rand();
    case "medium"
        scale = 0.25 + 0.45 * rand();
    case "boundary"
        scale = 0.70 + 0.25 * rand();
    case "quality_hard"
        scale = 0.35 + 0.45 * rand();
    otherwise
        scale = 0.25 + 0.45 * rand();
end

cap_p1 = 100;
cap_p2 = 200;

if mode == "quality_hard"
    C_pool = 1.42 + 0.07 * rand();
else
    C_pool = 1.05 + 0.35 * rand();
end
C_pool = min(max(C_pool, 1.0 + 1e-8), 1.5 - 1e-8);

total_p1 = cap_p1 * scale * (0.40 + 0.60 * rand());
total_p2 = cap_p2 * scale * (0.40 + 0.60 * rand());


split_p1 = 0.20 + 0.70 * rand();
f_o_p1 = split_p1 * total_p1;
f4_p1 = total_p1 - f_o_p1;


f_o_p2 = total_p2;
f4_p2 = 0;

pool_total = f_o_p1 + f_o_p2;


[f1_o, f2_o] = split_two_to_quality(pool_total, C_pool, 3.0, 1.0);
f3_o = 0;


x1 = f1_o;
x2 = f2_o;
x3 = f3_o;
x4 = f4_p1 + f4_p2;

x0 = [
    f1_o
    f2_o
    f3_o
    f_o_p1
    f_o_p2
    C_pool
    x1
    x2
    x3
    x4
    f4_p1
    f4_p2
    ];
end

function [high_flow, low_flow] = split_two_to_quality(total, q_target, q_high, q_low)
if total <= 1e-12
    high_flow = 0;
    low_flow = 0;
    return
end

alpha = (q_target - q_low) / (q_high - q_low);
alpha = min(max(alpha, 0), 1);

high_flow = alpha * total;
low_flow = total - high_flow;
end
