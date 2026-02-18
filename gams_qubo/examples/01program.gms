set i /1*5/;

parameter c(i)  /1 6,2 4,3 8,4 5,5 5/
          a1(i) /1 2,2 2,3 4,4 3,5 2/
          a2(i) /1 1,2 2,3 2,4 1,5 2/
          a3(i) /1 3,2 3,3 3,4 4,5 4/;

binary variable x(i);

free variable z;

equations obj, c1, c2, c3;

obj.. sum(i, c(i)*x(i)) =E= z;

c1.. sum(i, a1(i)*x(i)) =L= 7;

c2.. sum(i, a2(i)*x(i)) =E= 4;

c3.. sum(i, a3(i)*x(i)) =G= 5;

Model zeroOneProblem /all/;

*Solve zeroOneProblem use mip max z;

option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' zeroOneProblem MIP max z 10

display x.l, z.l;