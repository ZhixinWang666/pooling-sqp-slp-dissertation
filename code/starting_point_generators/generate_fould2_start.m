function x0 = generate_fould2_start(mode)

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

caps = [100; 200; 100; 200];
quality_margin = 1e-8;


if mode == "quality_hard"
    C1_pool1 = 1.42 + 0.07 * rand();  
    C1_pool2 = 1.85 + 0.10 * rand();  
else
    C1_pool1 = 1.10 + 0.35 * rand();
    C1_pool2 = 1.55 + 0.35 * rand();
end

C1_pool1 = min(max(C1_pool1, 1.0 + quality_margin), 1.5 - quality_margin);
C1_pool2 = min(max(C1_pool2, 1.5 + quality_margin), 2.0 - quality_margin);

o1 = zeros(4,1);
o2 = zeros(4,1);
f3 = zeros(4,1);
f6 = zeros(4,1);

for p = 1:4
    total_p = caps(p) * scale * (0.40 + 0.60 * rand());

    switch p
        case 1
            vals = split_random(total_p, 4);
            o1(p) = vals(1);
            o2(p) = vals(2);
            f3(p) = vals(3);
            f6(p) = vals(4);

        case 2
            o1(p) = total_p;

        case 3
            vals = split_random(total_p, 4);
            o1(p) = vals(1);
            o2(p) = vals(2);
            f3(p) = vals(3);
            f6(p) = vals(4);

        case 4
            vals = split_random(total_p, 3);
            o1(p) = vals(1);
            o2(p) = vals(2);
            f3(p) = vals(3);
    end
end

T1 = sum(o1);
T2 = sum(o2);

[f1_o1, f2_o1] = split_two_to_quality(T1, C1_pool1, 3.0, 1.0);
[f4_o2, f5_o2] = split_two_to_quality(T2, C1_pool2, 3.5, 1.5);

x1 = f1_o1;
x2 = f2_o1;
x3 = sum(f3);
x4 = f4_o2;
x5 = f5_o2;
x6 = sum(f6);

x0 = [
    f1_o1
    f2_o1
    f4_o2
    f5_o2

    o1

    o2

    C1_pool1
    C1_pool2

    x1
    x2
    x3
    x4
    x5
    x6

    f3

    f6
    ];
end

function y = split_random(total, n)
if total <= 1e-12
    y = zeros(n,1);
    return
end

w = rand(n,1);
w = w / sum(w);
y = total * w;
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
