function x0 = generate_adhya4_start(mode)


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
    q_c1_o1 = 0.09 + 0.03 * rand();
else
    q_c1_o1 = 0.06 + 0.08 * rand();
end
q_c2_o1 = 0;
q_c3_o1 = 0;
q_c4_o1 = 1 - q_c1_o1;


if mode == "quality_hard"
    q_c6_o2 = 0.62 + 0.08 * rand();
    q_c7_o2 = 0;
else
    q_c6_o2 = 0.58 + 0.17 * rand();
    q_c7_o2 = 0.02 * rand();
end
q_c5_o2 = 0;
q_c8_o2 = 1 - q_c5_o2 - q_c6_o2 - q_c7_o2;
if q_c8_o2 < 0.20
    q_c8_o2 = 0.20;
    q_c6_o2 = 1 - q_c7_o2 - q_c8_o2;
end


cap_p2 = 25;
cap_p3 = 10;
cap_p5 = 15;


y_o1 = zeros(5,1);
y_o2 = zeros(5,1);

total_p2 = cap_p2 * scale * (0.50 + 0.50 * rand());
total_p3 = cap_p3 * scale * (0.50 + 0.50 * rand());
total_p5 = cap_p5 * scale * (0.30 + 0.70 * rand());

y_o1(2) = total_p2;
y_o2(3) = total_p3;

split_p5 = 0.30 + 0.40 * rand();
y_o1(5) = split_p5 * total_p5;
y_o2(5) = total_p5 - y_o1(5);


y_o1 = min(y_o1, [15; 25; 10; 20; 15]);
y_o2 = min(y_o2, [10; 20; 15; 15; 25]);


if sum(y_o1) > 85
    y_o1 = y_o1 * (85 / sum(y_o1));
end
if sum(y_o2) > 85
    y_o2 = y_o2 * (85 / sum(y_o2));
end

p = y_o1 + y_o2;

x0 = [
    q_c1_o1
    q_c2_o1
    q_c3_o1
    q_c4_o1
    q_c5_o2
    q_c6_o2
    q_c7_o2
    q_c8_o2

    y_o1

    y_o2

    p
    ];
assert(numel(x0) == 23, 'adhya4 start vector must have 23 variables.');
end
