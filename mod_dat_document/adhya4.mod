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

param rawcomp {I,Q} default 0;
param prodcomp_ub {K,Q};

var q {i in I, j in J} >= 0, <= poolconn[i,j]
    := if poolconn[i,j] > 0 then 0.1 else 0;
var y {j in J, k in K} >= 0, <= maxpoolprod[j,k] * poolprodconn[j,k]
    := if poolprodconn[j,k] > 0 then 0.1 else 0;
var p {k in K} >= 0
    := sum{j in J: poolprodconn[j,k] > 0} 0.1;

minimize Profit:
    sum{i in I} raw_price[i] *
        sum{j in J} q[i,j] * sum{k in K} y[j,k]
    - sum{k in K} prod_price[k] * p[k];

subject to PoolFraction {j in J}:
    sum{i in I} q[i,j] = 1;

subject to RawLimit {i in I}:
    raw_lo[i] <= sum{j in J} q[i,j] * sum{k in K} y[j,k] <= raw_up[i];

subject to PoolCap {j in J}:
    sum{k in K} y[j,k] <= poolcap[j];

subject to ProdAmt {k in K}:
    p[k] = sum{j in J} y[j,k];

subject to ProdBound {k in K}:
    prod_lo[k] <= p[k] <= prod_up[k];

subject to CompLimit {k in K, ql in Q}:
    sum{j in J}
        (sum{i in I} rawcomp[i,ql] * q[i,j]) * y[j,k]
    <= prodcomp_ub[k,ql] * p[k];
