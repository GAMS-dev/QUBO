set i /b1*b6/;

Parameter c(i) /b1 3,b2 2,b3 1,b4 1,b5 3,b6 2/;

binary variables x(i);

free variable z;

Equations obj, c1, c2, c3, c4;

obj.. sum(i, c(i)*x(i)) =E= z;

c1.. x('b1') + x('b3') + x('b6') =E= 1;

c2.. x('b2') + x('b3') + x('b5') + x('b6') =E= 1;

c3.. x('b3') + x('b4') + x('b5') =E= 1;

c4.. x('b1') + x('b2')+  x('b4') + x('b6') =E= 1;

Model setPartition /all/;

*Solve setPartition use MIP min z;

option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' setPartition MIP min z 10

display x.l, z.l;
