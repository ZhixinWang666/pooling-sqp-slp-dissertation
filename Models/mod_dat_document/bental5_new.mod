var x1 >= 0;         
var x2 >= 0;         
var x3 >= 0;          
var x4 >= 0;          
var x5 >= 0;         
var x6 >= 0;        
var x7 >= 0;         
var x8 >= 0;         
var x9 >= 0;         
var x10 >= 0;         
var x11 >= 0;        
var x12 >= 0;         
var x13 >= 0;          

var f1_o1 >= 0;       
var f2_o1 >= 0;        
var f3_o1 >= 0;        
var f4_o1 >= 0;       
var f5_o2 >= 0;       
var f6_o2 >= 0;       
var f7_o2 >= 0;       
var f8_o2 >= 0;       
var f9_o3 >= 0;      
var f10_o3 >= 0;       
var f11_o3 >= 0;      
var f12_o3 >= 0;


var o1_p1 >= 0;
var o1_p2 >= 0;
var o1_p3 >= 0;
var o1_p4 >= 0;
var o1_p5 >= 0;
var o2_p1 >= 0;
var o2_p2 >= 0;
var o2_p3 >= 0;
var o2_p4 >= 0;
var o2_p5 >= 0;
var o3_p1 >= 0;
var o3_p2 >= 0;
var o3_p3 >= 0;
var o3_p4 >= 0;
var o3_p5 >= 0;



var f13_p1 >= 0;
var f13_p2 >= 0;
var f13_p3 >= 0;
var f13_p4 >= 0;
var f13_p5 >= 0;


var C1_pool1 >= 1, <= 3;
var C1_pool2 >= 1, <= 3;
var C1_pool3 >= 1, <= 3;
var C2_pool1 >= 1, <= 5;
var C2_pool2 >= 1, <= 5;
var C2_pool3 >= 1, <= 5;


minimize Profit:
    -18 * (o1_p1 + o2_p1 + o3_p1 + f13_p1) 
    -15 * (o1_p2 + o2_p2 + o3_p2 + f13_p2)
    -19 * (o1_p3 + o2_p3 + o3_p3 + f13_p3)
    -16 * (o1_p4 + o2_p4 + o3_p4 + f13_p4)
    -14 * (o1_p5 + o2_p5 + o3_p5 + f13_p5)
    +(6 * x1 + 16 * x2 + 15 * x3 + 12 * x4
    +6 * x5 + 16 * x6 + 15 * x7 + 12 * x8
    +6 * x9 + 16 * x10 + 15 * x11 + 12 * x12 + 10*x13);















s.t. Pool1_balance:
    f1_o1 + f2_o1 + f3_o1 + f4_o1 = o1_p1 + o1_p2 + o1_p3 + o1_p4 + o1_p5;
s.t. Pool2_balance:
    f5_o2 + f6_o2 + f7_o2 + f8_o2 = o2_p1 + o2_p2 + o2_p3 + o2_p4 + o2_p5;
s.t. Pool3_balance:
    f9_o3 + f10_o3 + f11_o3 + f12_o3 = o3_p1 + o3_p2 + o3_p3 + o3_p4 + o3_p5;


s.t. LB_x1: x1 <= 600;
s.t. LB_x2: x2 <= 600;
s.t. LB_x3: x3 <= 600;
s.t. LB_x4: x4 <= 600;
s.t. LB_x5: x5 <= 600;
s.t. LB_x6: x6 <= 600;
s.t. LB_x7: x7 <= 600;
s.t. LB_x8: x8 <= 600;
s.t. LB_x9: x9 <= 600;
s.t. LB_x10: x10 <= 600;
s.t. LB_x11: x11 <= 600;
s.t. LB_x12: x12 <= 600;
s.t. LB_x13: x13 <= 600;

s.t. Bal_x1: x1 = f1_o1;
s.t. Bal_x2: x2 = f2_o1;
s.t. Bal_x3: x3 = f3_o1;
s.t. Bal_x4: x4 = f4_o1;
s.t. Bal_x5: x5 = f5_o2;
s.t. Bal_x6: x6 = f6_o2;
s.t. Bal_x7: x7 = f7_o2;
s.t. Bal_x8: x8 = f8_o2;
s.t. Bal_x9: x9 = f9_o3;
s.t. Bal_x10: x10 = f10_o3;
s.t. Bal_x11: x11 = f11_o3;
s.t. Bal_x12: x12 = f12_o3;

s.t. Bal_x13: x13 = f13_p1 + f13_p2 + f13_p3 + f13_p4 + f13_p5;

s.t. Pool1Cap: 
     f1_o1 + f2_o1 + f3_o1 + f4_o1 <= 600;
s.t. Pool2Cap: 
     f5_o2 + f6_o2 + f7_o2 + f8_o2 <= 600;
s.t. Pool3Cap: 
     f9_o3 + f10_o3 + f11_o3 + f12_o3 <= 600;

s.t. Cap_o1_p1: o1_p1 <= 100;      
s.t. Cap_o1_p2: o1_p2 <= 200;
s.t. Cap_o1_p3: o1_p3 <= 100;
s.t. Cap_o1_p4: o1_p4 <= 100;
s.t. Cap_o1_p5: o1_p5 <= 100;

s.t. Cap_o2_p1: o2_p1 <= 100;
s.t. Cap_o2_p2: o2_p2 <= 100;      
s.t. Cap_o2_p3: o2_p3 <= 100;
s.t. Cap_o2_p4: o2_p4 <= 200;
s.t. Cap_o2_p5: o2_p5 <= 100;

s.t. Cap_o3_p1: o3_p1 <= 200;
s.t. Cap_o3_p2: o3_p2 <= 100;
s.t. Cap_o3_p3: o3_p3 <= 100;
s.t. Cap_o3_p4: o3_p4 <= 100;
s.t. Cap_o3_p5: o3_p5 <= 100;



s.t. Cap_f13_p1: f13_p1 <= 100;
s.t. Cap_f13_p2: f13_p2 <= 200;
s.t. Cap_f13_p3: f13_p3 <= 100;
s.t. Cap_f13_p4: f13_p4 <= 100;
s.t. Cap_f13_p5: f13_p5 <= 100;


s.t. Cap_y_p1: o1_p1 + o2_p1 + o3_p1 + f13_p1  <= 100;
s.t. Cap_y_p2: o1_p2 + o2_p2 + o3_p2 + f13_p2  <= 200;
s.t. Cap_y_p3: o1_p3 + o2_p3 + o3_p3 + f13_p3  <= 100;
s.t. Cap_y_p4: o1_p4 + o2_p4 + o3_p4 + f13_p4  <= 100;
s.t. Cap_y_p5: o1_p5 + o2_p5 + o3_p5 + f13_p5  <= 100;




s.t. Def_C1_pool1:
    C1_pool1 * (f1_o1 + f2_o1 + f3_o1 + f4_o1) = 3*f1_o1 + 1*f2_o1 + 1.2*f3_o1 + 1.5*f4_o1;
s.t. Def_C1_pool2:
    C1_pool2 * (f5_o2 + f6_o2 + f7_o2 + f8_o2 ) = 3*f5_o2 + 1*f6_o2 + 1.2*f7_o2 + 1.5*f8_o2;
s.t. Def_C1_pool3:
    C1_pool3 * (f9_o3 + f10_o3 + f11_o3 + f12_o3) = 3*f9_o3 + 1*f10_o3 + 1.2*f11_o3 + 1.5*f12_o3;
s.t. Def_C2_pool1:    
    C2_pool1 * (f1_o1 + f2_o1 + f3_o1 + f4_o1) = 1*f1_o1 + 3*f2_o1 + 5*f3_o1 + 2.5*f4_o1;
s.t. Def_C2_pool2:
    C2_pool2 * (f5_o2 + f6_o2 + f7_o2 + f8_o2) = 1*f5_o2 + 3*f6_o2 + 5*f7_o2 + 2.5*f8_o2;
s.t. Def_C2_pool3:
    C2_pool3 * (f9_o3 + f10_o3 + f11_o3 + f12_o3) = 1*f9_o3 + 3*f10_o3 + 5*f11_o3 + 2.5*f12_o3;

s.t. C1_Quality_p1:
    (2 * f13_p1 + C1_pool1 * o1_p1 + C1_pool2 * o2_p1 + C1_pool3 * o3_p1 ) <= 2.5 * (o1_p1 + o2_p1 + o3_p1 + f13_p1);
s.t. C2_Quality_p1:
    (2.5 * f13_p1 + C2_pool1 * o1_p1 + C2_pool2 * o2_p1 + C2_pool3 * o3_p1 ) <= 2 * (o1_p1 + o2_p1 + o3_p1 + f13_p1);

s.t. C1_Quality_p2:
    (2 * f13_p2 + C1_pool1 * o1_p2 + C1_pool2 * o2_p2 + C1_pool3 * o3_p2 ) <= 1.5 * (o1_p2 + o2_p2 + o3_p2 + f13_p2);
s.t. C2_Quality_p2:
    (2.5 * f13_p2 + C2_pool1 * o1_p2 + C2_pool2 * o2_p2 + C2_pool3 * o3_p2 ) <= 2.5 * (o1_p2 + o2_p2 + o3_p2 + f13_p2);

s.t. C1_Quality_p3:
    (2 * f13_p3 + C1_pool1 * o1_p3 + C1_pool2 * o2_p3 + C1_pool3 * o3_p3 ) <= 2 * (o1_p3 + o2_p3 + o3_p3 + f13_p3);
s.t. C2_Quality_p3:
    (2.5 * f13_p3 + C2_pool1 * o1_p3 + C2_pool2 * o2_p3 + C2_pool3 * o3_p3 ) <= 2.6 * (o1_p3 + o2_p3 + o3_p3 + f13_p3);
    
   
s.t. C1_Quality_p4:
    (2 * f13_p4 + C1_pool1 * o1_p4 + C1_pool2 * o2_p4 + C1_pool3 * o3_p4 ) <= 2 * (o1_p4 + o2_p4 + o3_p4 + f13_p4);
s.t. C2_Quality_p4:
    (2.5 * f13_p4 + C2_pool1 * o1_p4 + C2_pool2 * o2_p4 + C2_pool3 * o3_p4 ) <= 2 * (o1_p4 + o2_p4 + o3_p4 + f13_p4);

s.t. C1_Quality_p5:
    (2 * f13_p5 + C1_pool1 * o1_p5 + C1_pool2 * o2_p5 + C1_pool3 * o3_p5 ) <= 2 * (o1_p5 + o2_p5 + o3_p5 + f13_p5);
s.t. C2_Quality_p5:
    (2.5 * f13_p5 + C2_pool1 * o1_p5 + C2_pool2 * o2_p5 + C2_pool3 * o3_p5 ) <= 2 * (o1_p5 + o2_p5 + o3_p5 + f13_p5);

    
    
    
    
    
