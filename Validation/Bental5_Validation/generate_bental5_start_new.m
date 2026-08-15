function x0 = generate_bental5_start_new(mode)


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


if mode == "quality_hard"
    a1 = 0.50 + 0.03 * rand();
else
    a1 = 0.52 + 0.18 * rand();
end
a1 = min(max(a1, 0.50), 0.72);
q_o1 = [a1; 1 - a1; 0; 0];


if mode == "quality_hard"
    a2 = 0.005 + 0.020 * rand();
else
    a2 = 0.010 + 0.040 * rand();
end
a2 = min(a2, 0.20);
q_o2 = [a2; 3 * a2; 0; 1 - 4 * a2];


q_o3 = [0.5; 0.5; 0; 0];


prod_cap = [100; 200; 100; 100; 100];
y_o1_cap = [100; 200; 100; 100; 100];
y_o2_cap = [100; 100; 100; 200; 100];
y_o3_cap = [200; 100; 100; 100; 100];
pool_cap = 600;

y_o1 = zeros(5,1);
y_o2 = zeros(5,1);
y_o3 = zeros(5,1);


y_o1(1) = prod_cap(1) * scale * (0.50 + 0.50 * rand());
y_o2(2) = min(y_o2_cap(2), prod_cap(2)) * scale * (0.50 + 0.50 * rand());
y_o2(3) = prod_cap(3) * scale * (0.30 + 0.70 * rand());
y_o3(4) = prod_cap(4) * scale * (0.30 + 0.70 * rand());
y_o3(5) = prod_cap(5) * scale * (0.30 + 0.70 * rand());


z_c13 = zeros(5,1);

y_o1 = min(y_o1, y_o1_cap);
y_o2 = min(y_o2, y_o2_cap);
y_o3 = min(y_o3, y_o3_cap);

if sum(y_o1) > pool_cap
    y_o1 = y_o1 * (pool_cap / sum(y_o1));
end
if sum(y_o2) > pool_cap
    y_o2 = y_o2 * (pool_cap / sum(y_o2));
end
if sum(y_o3) > pool_cap
    y_o3 = y_o3 * (pool_cap / sum(y_o3));
end

p = y_o1 + y_o2 + y_o3 + z_c13;
p = min(p, prod_cap);

x0 = [
    q_o1
    q_o2
    q_o3

    y_o1
    y_o2
    y_o3

    z_c13

    p
    ];

assert(numel(x0) == 37, 'bental5_q start vector must have 37 variables.');
end
