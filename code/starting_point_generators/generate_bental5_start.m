function x0 = generate_bental5_start(mode)

if nargin < 1
    mode = "low";
end

switch mode
    case "low"
        scale = 0.01 + 0.09 * rand();
        bias = "balanced";
    case "medium"
        scale = 0.20 + 0.50 * rand();
        bias = "balanced";
    case "boundary"
        scale = 0.70 + 0.30 * rand();
        bias = "boundary";
    case "quality_hard"
        scale = 0.35 + 0.55 * rand();
        bias = "high_quality";
    otherwise
        scale = 0.20 + 0.50 * rand();
        bias = "balanced";
end

pool_cap = 600;
arc_caps = 100 * ones(5,1);


T1 = scale * pool_cap * rand();
T2 = scale * pool_cap * rand();
T3 = scale * pool_cap * rand();

f_o1 = split_four(T1, bias);  
f_o2 = split_four(T2, bias);  
f_o3 = split_four(T3, bias);  

[C1_pool1, C2_pool1] = pool_quality_bental5(f_o1);
[C1_pool2, C2_pool2] = pool_quality_bental5(f_o2);
[C1_pool3, C2_pool3] = pool_quality_bental5(f_o3);


C1_pool1 = min(max(C1_pool1, 1), 3);
C1_pool2 = min(max(C1_pool2, 1), 3);
C1_pool3 = min(max(C1_pool3, 1), 3);

C2_pool1 = min(max(C2_pool1, 1), 5);
C2_pool2 = min(max(C2_pool2, 1), 5);
C2_pool3 = min(max(C2_pool3, 1), 5);


o1 = split_products(T1, arc_caps);
o2 = split_products(T2, arc_caps);
o3 = split_products(T3, arc_caps);


f13 = arc_caps .* scale .* rand(5,1);


x1  = f_o1(1);
x2  = f_o1(2);
x3  = f_o1(3);
x4  = f_o1(4);

x5  = f_o2(1);
x6  = f_o2(2);
x7  = f_o2(3);
x8  = f_o2(4);

x9  = f_o3(1);
x10 = f_o3(2);
x11 = f_o3(3);
x12 = f_o3(4);

x13 = sum(f13);

x0 = [
    f_o1;
    f_o2;
    f_o3;

    o1;
    o2;
    o3;

    C1_pool1; C1_pool2; C1_pool3;
    C2_pool1; C2_pool2; C2_pool3;

    x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12; x13;

    f13
    ];

end

function f = split_four(T, bias)

if T <= 1e-12
    f = zeros(4,1);
    return
end

switch bias
    case "high_quality"
        w = [4; 0.5; 0.5; 0.5] .* rand(4,1);
    case "boundary"
        if rand() < 0.5
            w = 0.05 * rand(4,1);
            w(randi(4)) = 1;
        else
            w = rand(4,1);
            w(randi(4)) = 10;
        end
    otherwise
        w = rand(4,1);
end

if sum(w) <= 1e-12
    w = ones(4,1);
end

f = T * w / sum(w);
end

function [q1, q2] = pool_quality_bental5(f)

T = sum(f);

if T <= 1e-12
    q1 = 2;
    q2 = 2.5;
    return
end

q1_data = [3; 1; 1.2; 1.5];
q2_data = [1; 3; 5; 2.5];

q1 = dot(q1_data, f) / T;
q2 = dot(q2_data, f) / T;
end

function y = split_products(T, caps)

y = zeros(5,1);

if T <= 1e-12
    return
end

remaining = min(T, sum(caps));
free = true(5,1);
weights = rand(5,1);

while remaining > 1e-10 && any(free)
    idx = find(free);
    w = weights(idx) / sum(weights(idx));
    proposal = remaining * w;

    over = proposal > caps(idx);

    if ~any(over)
        y(idx) = y(idx) + proposal;
        break
    end

    over_idx = idx(over);
    y(over_idx) = caps(over_idx);
    remaining = remaining - sum(caps(over_idx));
    free(over_idx) = false;
end
end
