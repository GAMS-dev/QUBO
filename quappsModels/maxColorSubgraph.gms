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

*option miqcp=cplex;
*Solve mcs min TOTCOST using miqcp;


$set method classic
$set solver cplex
$set max_iter 1
$set timeLimit 60
$set num_threads 1
$set log_on 0

$batInclude qubo_solve mcs miqcp min TOTCOST 10 %method% %solver% %max_iter% %timeLimit% %num_threads% %log_on%

display X.l, TOTCOST.l;