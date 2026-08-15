# === Variables ===
var x1 >= 0;         
var x2 >= 0;         
var x3 >= 0;          
var x4 >= 0;          
var x5 >= 0;         
var x6 >= 0;  

var f1_o1 >= 0;       
var f2_o1 >= 0;  
       
       
var f4_o2 >= 0;       
var f5_o2 >= 0; 


var o1_p1 >= 0;    
var o1_p2 >= 0;     
var o1_p3 >= 0;     
var o1_p4 >= 0;    
var o2_p1 >= 0;    
var o2_p2 >= 0;     
var o2_p3 >= 0;     
var o2_p4 >= 0;      


var f3_p1 >= 0;      
var f3_p2 >= 0;      
var f3_p3 >= 0;       
var f3_p4 >= 0;       
var f6_p1 >= 0;      
var f6_p2 >= 0;      
var f6_p3 >= 0;       
var f6_p4 >= 0;       


var C1_pool1 >= 0, <= 3;
var C1_pool2 >= 0, <= 3.5;


minimize Profit:
    -9 * (o1_p1 + o2_p1 + f3_p1 + f6_p1) 
    -15 * (o1_p2 + o2_p2 + f3_p2 + f6_p2)
    -6 * (o1_p3 + o2_p3 + f3_p3 + f6_p3)
    -12 * (o1_p4 + o2_p4 + f3_p4 + f6_p4)
    +(6 * x1 + 16 * x2 + 10 * x3 + 3 * x4 + 13 * x5 + 7 * x6);
    
s.t. Pool1_balance:
    f1_o1 + f2_o1 = o1_p1 + o1_p2 + o1_p3 + o1_p4;
s.t. Pool2_balance:
    f4_o2 + f5_o2 = o2_p1 + o2_p2 + o2_p3 + o2_p4;


s.t. LB_x1: x1 <= 600;
s.t. LB_x2: x2 <= 600;
s.t. LB_x3: x3 <= 600;
s.t. LB_x4: x4 <= 600;
s.t. LB_x5: x5 <= 600;
s.t. LB_x6: x6 <= 600;



s.t. Bal_x1: x1 = f1_o1;
s.t. Bal_x2: x2 = f2_o1;
s.t. Bal_x4: x4 = f4_o2;
s.t. Bal_x5: x5 = f5_o2;

s.t. Bal_x3: x3 = f3_p1 + f3_p2 + f3_p3 + f3_p4;
s.t. Bal_x6: x6 = f6_p1 + f6_p2 + f6_p3 + f6_p4;



s.t. Pool1Cap: 
     f1_o1 + f2_o1  <= 600;
s.t. Pool2Cap: 
     f4_o2 + f5_o2  <= 600;



s.t. Cap_o1_p1: o1_p1 <= 100;      
s.t. Cap_o1_p2: o1_p2 <= 200;
s.t. Cap_o1_p3: o1_p3 <= 100;
s.t. Cap_o1_p4: o1_p4 <= 200;

s.t. Cap_o2_p1: o2_p1 <= 100;
s.t. Cap_o2_p2: o2_p2 <= 200;      
s.t. Cap_o2_p3: o2_p3 <= 100;
s.t. Cap_o2_p4: o2_p4 <= 200;


s.t. Cap_f3_p1: f3_p1 <= 100;
s.t. Cap_f3_p2: f3_p2 <= 200;
s.t. Cap_f3_p3: f3_p3 <= 100;
s.t. Cap_f3_p4: f3_p4 <= 200;

s.t. Cap_f6_p1: f6_p1 <= 100;
s.t. Cap_f6_p2: f6_p2 <= 200;
s.t. Cap_f6_p3: f6_p3 <= 100;
s.t. Cap_f6_p4: f6_p4 <= 200;



s.t. Cap_y_p1: o1_p1 + o2_p1 + f3_p1 + f6_p1  <= 100;
s.t. Cap_y_p2: o1_p2 + o2_p2 + f3_p2 + f6_p2  <= 200;
s.t. Cap_y_p3: o1_p3 + o2_p3 + f3_p3 + f6_p3  <= 100;
s.t. Cap_y_p4: o1_p4 + o2_p4 + f3_p4 + f6_p4  <= 200;


s.t. Def_C1_pool1:
    C1_pool1 * (f1_o1 + f2_o1) = 3*f1_o1 + 1*f2_o1;
s.t. Def_C1_pool2:
    C1_pool2 * (f4_o2 + f5_o2) = 3.5*f4_o2 + 1.5*f5_o2;
    

s.t. C1_Quality_p1:
    (2 * f3_p1 + 2.5 * f6_p1 +    C1_pool1 * o1_p1 + C1_pool2 * o2_p1 ) <= 2.5 * ( f3_p1 + f6_p1 + o1_p1 + o2_p1);


s.t. C1_Quality_p2:
    (2 * f3_p2 + 2.5 * f6_p2 +    C1_pool1 * o1_p2 + C1_pool2 * o2_p2 ) <= 1.5 * ( f3_p2 + f6_p2 + o1_p2 + o2_p2);


s.t. C1_Quality_p3:
    (2 * f3_p3 + 2.5 * f6_p3 +    C1_pool1 * o1_p3 + C1_pool2 * o2_p3 ) <= 3 * ( f3_p3 + f6_p3 + o1_p3 + o2_p3);
   
s.t. C1_Quality_p4:
    (2 * f3_p4 + 2.5 * f6_p4 +    C1_pool1 * o1_p4 + C1_pool2 * o2_p4 ) <= 2 * ( f3_p4 + f6_p4 + o1_p4 + o2_p4);
    






