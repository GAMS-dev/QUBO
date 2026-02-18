set i /1*4/;

Binary variable x(i);

Free variable z;

Equations obj, c1 , c2, c3;

obj.. sum(i, x(i)) =E= z;

c1.. sum(i$(not sameas(i,'2')), x(i)) =L= 1;

c2.. sum(i$(ord(i) <= 2), x(i)) =L= 1;

c3..  x('2') + x('3') =E= 1;

Model setPacking /all/;

*Solve setPacking using mip max z;

option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' setPacking MIP max z 6 -solver=cplex -timeLimit=60 -numThreads=2 -logOn=2

display x.l, z.l;