function x0 = generate_adhya3_start(mode)


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

if mode == "quality_hard"
    q_c6_o3 = 0.03 + 0.05 * rand();
    q_c7_o3 = 0.02 + 0.04 * rand();
else
    q_c6_o3 = 0.03 + 0.12 * rand();
    q_c7_o3 = 0.02 + 0.10 * rand();
end
if q_c6_o3 + q_c7_o3 > 0.30
    factor = 0.30 / (q_c6_o3 + q_c7_o3);
    q_c6_o3 = factor * q_c6_o3;
    q_c7_o3 = factor * q_c7_o3;
end
q_c8_o3 = 1 - q_c6_o3 - q_c7_o3;


cap_p2 = 25;
cap_p4 = 10;


y_o1 = zeros(4,1);
y_o2 = zeros(4,1);
y_o3 = zeros(4,1);

total_p2 = cap_p2 * scale * (0.50 + 0.50 * rand());
total_p4 = cap_p4 * scale * (0.50 + 0.50 * rand());


if mode == "quality_hard"
    o2_share_p2 = 0.02 * rand();
    o3_share_p2 = 0.20 + 0.15 * rand();
    o3_share_p4 = 0.20 + 0.20 * rand();
else
    o2_share_p2 = 0.05 * rand();
    o3_share_p2 = 0.15 + 0.25 * rand();
    o3_share_p4 = 0.15 + 0.30 * rand();
end
o1_share_p2 = 1 - o2_share_p2 - o3_share_p2;

o2_share_p4 = 0.10 + 0.20 * rand();
o1_share_p4 = 1 - o2_share_p4 - o3_share_p4;
if o1_share_p4 < 0.20
    o1_share_p4 = 0.20;
    remaining = 1 - o1_share_p4;
    total_other = o2_share_p4 + o3_share_p4;
    o2_share_p4 = remaining * o2_share_p4 / total_other;
    o3_share_p4 = remaining * o3_share_p4 / total_other;
end

y_o1(2) = o1_share_p2 * total_p2;
y_o2(2) = o2_share_p2 * total_p2;
y_o3(2) = o3_share_p2 * total_p2;

y_o1(4) = o1_share_p4 * total_p4;
y_o2(4) = o2_share_p4 * total_p4;
y_o3(4) = o3_share_p4 * total_p4;

y_o1 = min(y_o1, [10; 25; 30; 10]);
y_o2 = min(y_o2, [10; 10; 25; 30]);
y_o3 = min(y_o3, [30; 10; 10; 25]);

if sum(y_o1) > 75
    y_o1 = y_o1 * (75 / sum(y_o1));
end
if sum(y_o2) > 75
    y_o2 = y_o2 * (75 / sum(y_o2));
end
if sum(y_o3) > 75
    y_o3 = y_o3 * (75 / sum(y_o3));
end

p = y_o1 + y_o2 + y_o3;

x0 = [
    q_c1_o1
    q_c2_o1
    q_c3_o2
    q_c4_o2
    q_c5_o2
    q_c6_o3
    q_c7_o3
    q_c8_o3

    y_o1

    y_o2

    y_o3

    p
    ];
 assert(numel(x0) == 24, 'adhya3 start vector must have 24 variables.');
end
