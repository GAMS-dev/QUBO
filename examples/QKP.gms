set i /1*4/;

alias(i,j);

Parameter uc(i,j) /1.1 2,2.2 5, 3.3 2, 4.4 4,
                  1.2 8, 1.3 6, 1.4 10,
                  2.3 2, 2.4 6, 3.4 4/
          pc(i)  /1 8, 2 6, 3 5, 4 3/;
           

binary variable x(i);

free variable z;

equation obj, c1;

obj.. sum((i,j), x(i)*uc(i,j)*x(j)) =E= z;

c1.. sum(i, pc(i)*x(i)) =L= 16;

Model qkp /all/;

*option miqcp=cplex;
*solve qkp using miqcp max z;

option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' qkp miqcp max z 10

display x.l, z.l;