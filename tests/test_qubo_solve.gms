$if not set TESTDWAVE $set TESTDWAVE no

$log Test 1: No Continuous variables allowed
$onEcho>test1.gms
Positive Variable x1;
Free Variable z;
Equation obj, c1;
obj.. x1 =E= z;
c1.. x1 =E= 1; 
Model test1 /all/;
$batinclude %QUBO_PATH% test1 MIP max z 1
$offEcho

$call gams test1.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test1.gms
$call rm test1.lst
$call rm test1.gdx

$exit

$log Test 2: No Real-valued Coefficients allowed
$onEcho>test2.gms
Binary Variable x1;
Free Variable z;
Equation obj, c1;
obj.. x1 =E= z;
c1.. 2.3*x1 =E= 1;
Model test2 /all/;
$batinclude %QUBO_PATH% test2 MIP max z 1
$offEcho
$call gams test2.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test2.gms
$call rm test2.lst
$call rm test2.gdx

$log Test 3: Objective Function must be defined using a scalar equation
$onEcho>test3.gms
Set i/1*5/;
Binary Variable x(i);
Free variable z;
Equation obj(i), c1;
obj(i).. x(i) =G= z;
c1.. sum(i, x(i)) =E= 1; 
Model test3 /all/;
$batinclude %QUBO_PATH% test3 MIP max z 1
$offEcho
$call gams test3.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test3.gms
$call rm test3.lst
$call rm test3.gdx

$log Test 4(i): Infeasible constraint raises Exception
$onEcho>test4_1.gms
Set i/1*5/;
Binary Variable x(i);
Free variable z;
Equation obj, c1;
obj.. sum(i, x(i)) =E= z;
c1.. sum(i, x(i)) =E= 10; 
Model test4_1 /all/;
$batinclude %QUBO_PATH% test4_1 MIP max z 1
$offEcho
$call gams test4_1.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test4_1.gms
$call rm test4_1.lst
$call rm test4_1.gdx

$log Test 4(ii): Infeasible constraint raises Exception
$onEcho>test4_2.gms
Set i/1*5/;
Binary Variable x(i);
Free variable z;
Equation obj, c1;
obj.. sum(i, x(i)) =E= z;
c1.. sum(i, x(i)) =L= -10; 
Model test4_2 /all/;
$batinclude %QUBO_PATH% test4_2 MIP max z 1
$offEcho
$call gams test4_2.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test4_2.gms
$call rm test4_2.lst
$call rm test4_2.gdx

$log Test 4(iii): Infeasible constraint raises Exception
$onEcho>test4_3.gms
Set i/1*5/;
Binary Variable x(i);
Free variable z;
Equation obj, c1;
obj.. sum(i, x(i)) =E= z;
c1.. sum(i, x(i)) =G= 10;
Model test4_3 /all/;
$batinclude %QUBO_PATH% test4_3 MIP max z 1
$offEcho
$call gams test4_3.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test4_3.gms
$call rm test4_3.lst
$call rm test4_3.gdx

$log Test 5: Upper bound of Integer variables must be less than 1e4
$onEcho>test5.gms
Set i/1*5/;
Integer Variable x(i);
x.up(i) = 1e5;
Free variable z;
Equation obj, c1;
obj.. sum(i, x(i)) =E= z;
c1.. sum(i, x(i)) =L= 50;
Model test5 /all/;
$batinclude %QUBO_PATH% test5 MIP max z 1
$offEcho
$call gams test5.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test5.gms
$call rm test5.lst
$call rm test5.gdx

$log Test 6: Non-Linear constraints
$onEcho>test6.gms
Set i/1*5/;
alias (i,j);
Binary Variable x(i);
Free variable z;
Equation obj, c1;
obj.. sum(i, x(i)) =E= z;
c1.. sum((i,j), x(i)*x(j)) =G= 2;
Model test6 /all/;
$batinclude %QUBO_PATH% test6 MIQCP max z 1
$offEcho
$call gams test6.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$ifE errorLevel=0 $abort 'Expect error'
$call rm test6.gms
$call rm test6.lst
$call rm test6.gdx

$log Test 7: Check correctness of Reformulation
$onEcho>test7.gms
Set i/1*5/;
Binary Variable x(i);
parameter cost(i);
cost(i) = UniformInt(1,10);
Free variable z;
Equation obj, c1;
obj.. sum(i, cost(i)*x(i)) =E= z;
c1.. sum(i, x(i)) =L= 3;
Model test7 /all/;
$batinclude %QUBO_PATH% test7 MIQCP max z 5
if((sum(i, round(x.l(i))) ne 3) or (round(z.l) ne 19), abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test7.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test7.gms
$call rm test7.lst
$call rm test7.gdx
$call rm qout_test7.gdx

$log Test 8: Support for non indexed/Scalar GAMS file
$onEcho>test8.gms
Binary Variable x1, x2;
Free Variable z;
Equation obj, c1;
obj.. 2*x1 + 3*x2 =E= z;
c1.. x1+x2 =E= 1;
Model test8 /all/;
$batinclude %QUBO_PATH% test8 MIP max z 5
if(round(x2.l) ne 1 or round(z.l) ne 3, abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test8.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test8.gms
$call rm test8.lst
$call rm test8.gdx
$call rm qout_test8.gdx

$ifThenI.testqpu %TESTDWAVE%==yes
$log Test 9: Solving on Dwave QPU
$onEcho>test9.gms
Binary Variable x1, x2;
Free Variable z;
Equation obj, c1;
obj.. 2*x1 + 3*x2 =E= z;
c1.. x1+x2 =E= 1;
Model test9 /all/;
$batinclude %QUBO_PATH% test9 MIP max z 5 qpu "" 1 5
if(round(x2.l) ne 1 or round(z.l) ne 3, abort 'Solution is not right');
$offEcho
$call.checkErrorLevel gams test9.gms --QUBO_PATH %system.fp%..%system.dirsep%qubo_solve.gms
$call rm test9.gms
$call rm test9.lst
$call rm test9.gdx
$endif.testqpu




$libInclude moo.gms 