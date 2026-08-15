# === Sets ===
set C;                       
set P;                        
set O;                        
set CP_LINK within {C,O};     

# === Parameters ===
param ub_c{C};                
param cost_c{C};              
param conc_c{C};              

param ub_p{P};               
param price_p{P};             
param target_conc_p{P};      

param cap_o{O};             
param cap_op{P,O};            

# === Variables ===
var x{C} >= 0;                        
var f{(c,o) in CP_LINK} >= 0;         
var o2p{O,P} >= 0;                    
var C1_o{O} >= 0, <= max{c in C} conc_c[c];  

# === Objective===
minimize NegProfit:
    -sum{p in P} price_p[p] * sum{o in O} o2p[o,p]
   + sum{c in C} cost_c[c] * x[c];


s.t. Bal_x{c in C}:
    x[c] = sum{o in O: (c,o) in CP_LINK} f[c,o];


s.t. PoolConc{o in O}:
    C1_o[o] * sum{c in C: (c,o) in CP_LINK} f[c,o]
  = sum{c in C: (c,o) in CP_LINK} conc_c[c] * f[c,o];


s.t. PoolBalance{o in O}:
    sum{c in C: (c,o) in CP_LINK} f[c,o] = sum{p in P} o2p[o,p];


s.t. PoolCap{o in O}:
    sum{c in C: (c,o) in CP_LINK} f[c,o] <= cap_o[o];


s.t. Cap_op{o in O, p in P}:
    o2p[o,p] <= cap_op[p,o];   


s.t. ProdCap{p in P}:
    sum{o in O} o2p[o,p] <= ub_p[p];


s.t. ProdConc{p in P}:
    sum{o in O} C1_o[o] * o2p[o,p] <= target_conc_p[p] * sum{o in O} o2p[o,p];
