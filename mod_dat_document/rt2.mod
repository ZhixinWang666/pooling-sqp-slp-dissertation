set I;  # raw materials
set J;  # pools
set K;  # products
set Q;  # qualities

param raw_lo {I};
param raw_up {I};
param raw_price {I};

param poolcap {J};
param poolconn {I,J};  

param prod_lo {K};
param prod_up {K};
param prod_price {K};

param poolprodconn {J,K};     
param maxpoolprod {J,K};      
param maxrawprod {I,K};      

param rawcomp {I,Q};
param prodcomp_lo {K,Q};
param prodcomp_ub {K,Q};

var q {i in I, j in J} >= 0, <= poolconn[i,j]
    := if poolconn[i,j] > 0 then 0.1 else 0;
var y {j in J, k in K} >= 0, <= maxpoolprod[j,k] * poolprodconn[j,k]
    := if poolprodconn[j,k] > 0 then 0.1 else 0;
var z {i in I, k in K} >= 0, <= maxrawprod[i,k]
    := if maxrawprod[i,k] > 0 then 0.1 else 0;
var p {k in K} >= 0
    := sum{j in J: poolprodconn[j,k] > 0} 0.1
     + sum{i in I: maxrawprod[i,k] > 0} 0.1;

minimize Profit:
    sum{i in I} raw_price[i] *
        (sum{k in K} z[i,k] + sum{j in J} q[i,j] * sum{k in K} y[j,k])
    - sum{k in K} prod_price[k] * p[k];

subject to PoolFraction {j in J}:
    sum{i in I} q[i,j] = 1;

subject to RawLimit {i in I}:
    raw_lo[i] <= sum{k in K} z[i,k]
               + sum{j in J} q[i,j] * sum{k in K} y[j,k]
               <= raw_up[i];

subject to PoolCap {j in J}:
    sum{k in K} y[j,k] <= poolcap[j];

subject to ProdAmt {k in K}:
    p[k] = sum{j in J} y[j,k] + sum{i in I} z[i,k];

subject to ProdBound {k in K}:
    prod_lo[k] <= p[k] <= prod_up[k];

subject to CompLower {k in K, ql in Q}:
    sum{j in J}
        (sum{i in I} rawcomp[i,ql] * q[i,j]) * y[j,k]
    + sum{i in I} rawcomp[i,ql] * z[i,k]
    >= prodcomp_lo[k,ql] * p[k];

subject to CompUpper {k in K, ql in Q}:
    sum{j in J}
        (sum{i in I} rawcomp[i,ql] * q[i,j]) * y[j,k]
    + sum{i in I} rawcomp[i,ql] * z[i,k]
    <= prodcomp_ub[k,ql] * p[k];
