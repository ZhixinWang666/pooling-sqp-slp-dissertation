function x0 = generate_adhya1_start(mode)


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
    q_c1_o1 = 0.835 + 0.020 * rand();
else
    q_c1_o1 = 0.84 + 0.05 * rand();
end
q_c1_o1 = min(max(q_c1_o1, 0.835), 0.895);
q_c2_o1 = 1 - q_c1_o1;


if mode == "quality_hard"
    q_c4_o2 = 0.15 + 0.20 * rand();
else
    q_c4_o2 = 0.05 + 0.30 * rand();
end
q_c3_o2 = 0;
q_c5_o2 = 1 - q_c4_o2;

cap_p3 = 30;
cap_p4 = 10;


y_o1 = zeros(4,1);
y_o2 = zeros(4,1);

y_o1(3) = cap_p3 * scale * (0.50 + 0.50 * rand());
y_o2(4) = cap_p4 * scale * (0.50 + 0.50 * rand());


pool1_total = sum(y_o1);
pool2_total = sum(y_o2);
if pool1_total > 75
    y_o1 = y_o1 * (75 / pool1_total);
end
if pool2_total > 75
    y_o2 = y_o2 * (75 / pool2_total);
end

p = y_o1 + y_o2;

x0 = [
    q_c1_o1
    q_c2_o1
    q_c3_o2
    q_c4_o2
    q_c5_o2

    y_o1

    y_o2

    p
    ];
assert(numel(x0) == 17, 'adhya1 start vector must have 17 variables.');
end
