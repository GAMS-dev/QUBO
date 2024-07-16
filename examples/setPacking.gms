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

$set method classic
$set solver cplex
$set max_iter 1
$set timeLimit 60
$set num_threads 3
$set log_on 2

$batInclude qubo_solve setPacking MIP max z 6 %method% %solver% %max_iter% %timeLimit% %num_threads% %log_on%

display x.l, z.l;