Set nodes /a,b,c,d/
    num_clr /0,1,2/;
    
alias (n,n1,n2,nodes), (c,c1,c2,num_clr);

Set edges(n1,n2) /a.b, b.c, b.d, c.d/;

Binary Variable X(n,c)  '1 if node n is colored with color c and 0 otherwise';

Free Variable TOTCOST;

Equations
    cost_fn             'Objective function',
    eq_get_one_clr(n)   'Each node gets only one color';
    
cost_fn..               TOTCOST =E= sum((edges(n1,n2), c), X(n1,c)*X(n2,c));

eq_get_one_clr(n)..     sum(c, X(n,c)) =E= 1;

Model mcs /all/;

option limrow=0, limcol=0;

*option miqcp=cplex;
*Solve mcs min TOTCOST using miqcp;

$batInclude '..\qubo_solve.gms' mcs miqcp min TOTCOST 10 -solver=cplex -timeLimit=60

display X.l, TOTCOST.l;