* A random IP to test the qubo reformulation in a more general setting*

set i /b1*b4/;

parameter cost(i) /b1 2, b2 4, b3 7, b4 9/;

integer Variables x(i);

loop(i,
    x.up(i) = 9;
);

Binary Variable y(i);

Free Variables z;

Integer Variables newX(i);

loop(i,
    newX.up(i) = 5;
);

Equations obj, c1, c2, c3;

obj.. sum(i, cost(i)*x(i)) + sum(i, cost(i)*y(i)) - sum(i, newX(i)) =E= z;

c1.. sum(i, x(i))       =L= 20;
c2.. sum(i, y(i))       =L= 3;
c3.. sum(i, newX(i))    =G= 3;

model demo_model /all/;

*solve demo_model using mip max z;
option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' demo_model mip max z 10