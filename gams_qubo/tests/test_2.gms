$if not set TESTDWAVE $set TESTDWAVE no

$log Test 1: Check correctness of Reformulation
$onEcho>test1.gms
Set i/1*5/;
Binary Variable x(i);
parameter cost(i);
cost(i) = UniformInt(1,10);
Free variable z;
Equation obj, c1;
obj.. sum(i, cost(i)*x(i)) =E= z;
c1.. sum(i, x(i)) =L= 3;
Model test1 /all/;
$batinclude %QUBO_PATH% test1 MIQCP max z 5
if((sum(i, round(x.l(i))) ne 3) or (round(z.l) ne 19), abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test1.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test1.gms
$call rm test1.lst
$call rm test1.gdx
$call rm qout_test1.gdx

$log Test 2: Support for non indexed/Scalar GAMS file
$onEcho>test2.gms
Binary Variable x1, x2;
Free Variable z;
Equation obj, c1;
obj.. 2*x1 + 3*x2 =E= z;
c1.. x1+x2 =E= 1;
Model test2 /all/;
$batinclude %QUBO_PATH% test2 MIP max z 5
if(round(x2.l) ne 1 or round(z.l) ne 3, abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test2.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test2.gms
$call rm test2.lst
$call rm test2.gdx
$call rm qout_test2.gdx

$ifThenI.testqpu %TESTDWAVE%==yes
$log Test 3: Solving on Dwave QPU
$onEcho>test3.gms
Binary Variable x1, x2;
Free Variable z;
Equation obj, c1;
obj.. 2*x1 + 3*x2 =E= z;
c1.. x1+x2 =E= 1;
Model test3 /all/;
$batinclude %QUBO_PATH% test3 MIP max z 5 -method=qpu -timeLimit=5
if(round(x2.l) ne 1 or round(z.l) ne 3, abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test3.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test3.gms
$call rm test3.lst
$call rm test3.gdx
$endif.testqpu