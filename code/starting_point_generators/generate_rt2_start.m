function x0 = generate_rt2_start(mode)


if nargin < 1
    mode = "low";
else
    mode = string(mode);
end

switch mode
    case "low"
        scale = 0.10 + 0.20 * rand();
    case "medium"
        scale = 0.35 + 0.35 * rand();
    case "boundary"
        scale = 0.65 + 0.25 * rand();
    case "quality_hard"
        scale = 0.20 + 0.25 * rand();
    otherwise
        scale = 0.25 + 0.35 * rand();
end


if mode == "quality_hard"
    q_c1_o1 = 0.70 + 0.03 * rand();
else
    q_c1_o1 = 0.68 + 0.07 * rand();
end
q_c3_o1 = 0;
q_c2_o1 = 1 - q_c1_o1 - q_c3_o1;


if mode == "quality_hard"
    q_c1_o2 = 0.285 + 0.010 * rand();
    q_c3_o2 = 0.490 + 0.035 * rand();
else
    q_c1_o2 = 0.270 + 0.025 * rand();
    q_c3_o2 = 0.470 + 0.060 * rand();
end
q_c2_o2 = 1 - q_c1_o2 - q_c3_o2;


if q_c2_o2 < 0.05
    q_c2_o2 = 0.05;
    q_c3_o2 = 1 - q_c1_o2 - q_c2_o2;
end

y_o1 = zeros(3, 1);
y_o2 = zeros(3, 1);


total_o1 = 10 + 2.2 * scale;
split_p1 = 0.46 + 0.08 * rand();
y_o1(1) = split_p1 * total_o1;
y_o1(3) = total_o1 - y_o1(1);

if y_o1(1) < 5
    y_o1(1) = 5;
    y_o1(3) = total_o1 - y_o1(1);
end
if y_o1(3) < 5
    y_o1(3) = 5;
    y_o1(1) = total_o1 - y_o1(3);
end


y_o2(2) = 5 + 1.2 * scale * (0.5 + 0.5 * rand());
y_o2(2) = min(y_o2(2), 0.95 * 5 / q_c3_o2);


z_c1_p2 = 0;
z_c2_p1 = 0;
z_c2_p3 = 0;
z_c3_p1 = 0;

p1 = y_o1(1) + y_o2(1) + z_c2_p1 + z_c3_p1;
p2 = y_o1(2) + y_o2(2) + z_c1_p2;
p3 = y_o1(3) + y_o2(3) + z_c2_p3;

x0 = [
    q_c1_o1
    q_c1_o2
    q_c2_o1
    q_c2_o2
    q_c3_o1
    q_c3_o2
    y_o1
    y_o2
    z_c1_p2
    z_c2_p1
    z_c2_p3
    z_c3_p1
    p1
    p2
    p3
];
assert(numel(x0) == 19, 'rt2 start vector must have 19 variables.');

end
