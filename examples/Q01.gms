$onText
This example handles case where there are quadratic constraints and objective function.

It should be noted that the quadratic constraints do not contain pair-wise quadratic terms. This is a limitation of QUBO reformulation. The quadratic terms otherwise, for e.g., (x_1)^2 can be treated as x_1 since x_1 is binary.
$offText

set i /1*5/;

alias (i,j);

parameter c(i)  /1 6,2 4,3 8,4 5,5 5/
          a1(i) /1 2,2 2,3 4,4 3,5 2/
          a2(i) /1 1,2 2,3 2,4 1,5 2/
          a3(i) /1 3,2 3,3 3,4 4,5 4/;

binary variable x(i);

free variable z;

equations obj, c1, c2, c3;

obj.. sum(i, x(i)*c(i)*x(i)) =E= z;

c1.. sum(i, a1(i)*x(i)) =L= 12;

c2.. sum(i$[ord(i) <=3], a2(i)*x(i)) - sum(i$[ord(i) > 3], a2(i)*x(i)) =G= 5;

c3.. sum(i$[ord(i)<= 3], a3(i)*x(i))=L= 10;

Model quadZeroOne /all/;

*option miqcp=cplex;
*Solve quadZeroOne use miqcp min z;

$batInclude qubo_solve quadZeroOne MIQCP min z 10

display x.l, z.l;