var x1 >= 0;
var x2 >= 0;
var x3 >= 0;
var x4 >= 0;

var f1_o >= 0;
var f2_o >= 0;
var f3_o >= 0;

var f_o_p1 >= 0;
var f_o_p2 >= 0;

var f4_p1 >= 0;
var f4_p2 >= 0;


minimize Profit:
    -9 * (f_o_p1 +f4_p1) - 15 * (f_o_p2 + f4_p2)
     +(6 * x1 + 15 * x2 + 16 * x3 + 10 * x4);

s.t. Pool_balance:
    f1_o + f2_o + f3_o = f_o_p1 + f_o_p2;


s.t. LB_x1: x1 <= 300;
s.t. LB_x2: x2 <= 50;
s.t. LB_x3: x3 <= 300;
s.t. LB_x4: x4 <= 300;

s.t. Bal_x1: x1 = f1_o;
s.t. Bal_x2: x2 = f2_o;
s.t. Bal_x3: x3 = f3_o;
s.t. Bal_x4: x4 = f4_p1 + f4_p2;

s.t. PoolCap: f1_o + f2_o + f3_o <= 300;

s.t. Flow_o_p1: f_o_p1 <= 100;
s.t. Flow_o_p2: f_o_p2 <= 200;

s.t. Flow_4_p1: f4_p1 <= 100;
s.t. Flow_4_p2: f4_p2 <= 200;


s.t. Cap_y_p1: f_o_p1 + f4_p1 <= 100;
s.t. Cap_y_p2: f_o_p2 + f4_p2 <= 200;

var C_pool >= 0, <= 3;
s.t. Def_C_pool:
    C_pool * (f1_o + f2_o + f3_o) = 3*f1_o + 1*f2_o + 1*f3_o;

s.t. Quality_p1:

    (2 * f4_p1 + C_pool * f_o_p1) <= 2.5 * (f_o_p1 + f4_p1);


s.t. Quality_p2:

    (2 * f4_p2 + C_pool * f_o_p2) <= 1.5 * (f_o_p2 + f4_p2);
