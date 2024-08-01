Sets
    Nodes       'nodes' /a,b,c,d/
    Position    /0,1,2,3,4/;

alias (n,n1,n2,Nodes), (i,i1,position);

Set Edges(n1,n2) 'possible edges';

Parameter weights(n1,n2) / a.b 8, a.c 1, a.d 1,
                           b.a 4, b.c 1, b.d 6,
                           c.a 2, c.b 4, c.d 9,
                           d.a 7, d.b 1, d.c 5/;

Edges(n1,n2) = weights(n1,n2);

Free Variable TOTCOST;

Binary Variable X(n,i);

X.fx(n,i)$(n.first and i.first)                    = 1;
X.fx(n,i)$(n.first and i.last)                     = 1;
X.fx(n,i)$(n.first and not i.first and not i.last) = 0;
X.fx(n,i)$(not n.first and i.first)                = 0;
X.fx(n,i)$(not n.first and i.last)                 = 0;

Equations
    cost_fn                 'Objective function'
    eq_exact_pos(i)         'each position is used exactly once',
    eq_exact_node(n)        'visit each node exactly once';
    

cost_fn.. TOTCOST =E= sum((edges(n1,n2), i)$[not i.last], weights(n1,n2)*X(n1,i)*X(n2,i+1));

eq_exact_pos(i)$[not i.first and not i.last]..          sum(n$(not n.first), X(n,i)) =E= 1;

eq_exact_node(n)$[not n.first]..                        sum(i$[not i.first and not i.last], X(n,i)) =E= 1;


Model tsp /all/;

*option miqcp=cplex;
*Solve tsp using miqcp min TOTCOST;

* model attribute holdfixed results in fixed variables being treated as constants
tsp.holdfixed = 1;

$batInclude '..\qubo_solve.gms' tsp miqcp min TOTCOST 10 -solver=cplex -timeLimit=60

display X.l;