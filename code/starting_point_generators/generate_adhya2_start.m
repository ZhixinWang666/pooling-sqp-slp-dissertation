function x0 = generate_adhya2_start(mode)


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
    q_c1_o1 = 0.26 + 0.04 * rand();
else
    q_c1_o1 = 0.24 + 0.08 * rand();
end
q_c1_o1 = min(max(q_c1_o1, 0.22), 0.34);
q_c2_o1 = 1 - q_c1_o1;


if mode == "quality_hard"
    q_c3_o2 = 0;
    q_c4_o2 = 0.60 + 0.08 * rand();
else
    q_c3_o2 = 0.02 * rand();
    q_c4_o2 = 0.55 + 0.18 * rand();
end
q_c5_o2 = 1 - q_c3_o2 - q_c4_o2;
if q_c5_o2 < 0.10
    q_c5_o2 = 0.10;
    q_c4_o2 = 1 - q_c3_o2 - q_c5_o2;
end


cap_p2 = 25;
cap_p4 = 10;

y_o1 = zeros(4,1);
y_o2 = zeros(4,1);

total_p2 = cap_p2 * scale * (0.50 + 0.50 * rand());
total_p4 = cap_p4 * scale * (0.50 + 0.50 * rand());


if mode == "quality_hard"
    split_p2 = 0.70 + 0.10 * rand();
    split_p4 = 0.75 + 0.10 * rand();
else
    split_p2 = 0.62 + 0.22 * rand();
    split_p4 = 0.65 + 0.25 * rand();
end

y_o1(2) = split_p2 * total_p2;
y_o2(2) = total_p2 - y_o1(2);
y_o1(4) = split_p4 * total_p4;
y_o2(4) = total_p4 - y_o1(4);


y_o1 = min(y_o1, [10; 25; 30; 10]);
y_o2 = min(y_o2, [30; 10; 10; 25]);


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
assert(numel(x0) == 17, 'adhya2 start vector must have 17 variables.');
end
