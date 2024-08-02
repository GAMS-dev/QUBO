$eolCom !!

Sets F 'set of flights' /F1,F2,F3/,
     G 'set of gates'   /G1,G2/;     

alias (i,j,F), (l,m,G);

Scalar buffer_time 'the buffer time between two flights at the same gate' /10/;

Parameter
    arr_time(G)             'the time it takes for a passenger to get from gate G to baggage claim'   /G1 2, G2 2.5/,
    dep_time(G)             'the time it takes for a passenger to get from check-in to gate G'        /G1 3, G2 2/,
    trnsfr_time(l,m)        'the time it takes to get from gate l to gate m'                          /G1.G2 7, G2.G1 9, G1.G1 5, G2.G2 5/,
    passengers_in(F)        'the number of passengers from flight F which arrive at the airport'      /F1 50, F3 100/,
    passengers_out(F)       'the number of passengers which depart from the airport on flight F'      /F1 10, F2 25, F3 65/,
    passenger_trnsfr(i,j)   'the number of lay over passengers which arrive on flight i and depart with flight j'
                                                        /F1.F1 10, F1.F3 20, F2.F1 20, F2.F2 10, F3.F1 20, F3.F3 10/,
    time_in(F)              'the arrival time of flight F'                              /F2 20, F3 35/,
    time_out(F)             'the departure time of flight F'                            /F1 16, F2 30, F3 50/;
    

Binary Variables X(F,G) '1 iff flight F is assigned to gate G, 0 otherwise';

Free Variable TOTCOST;

Equations
    cost_fn                     'Objective Function',
    eq_use_one_gate(i)          'Allow only one gate per flight',
    eq_restrict_arrival(i,j,l)  'NO arrival before departure at the same gate'
    eq_restrict_arrival_linear(i,j,l) 'linear alternative for its quadratic counterpart';
    
Set P(i,j);

P(i,j)$((time_in(i) < time_in(j)) and (time_in(j) < time_out(i) + buffer_time)) = yes;

cost_fn..                           TOTCOST =E= sum((i,l), (passengers_out(i)*dep_time(l) + passengers_in(i)*arr_time(l))*X(i,l)) +
                                                sum((i,j,l,m), passenger_trnsfr(i,j)*trnsfr_time(l,m)*X(i,l)*X(j,m));

eq_use_one_gate(i)..                sum(l, X(i,l)) =E= 1;

eq_restrict_arrival(i,j,l)$(P(i,j)).. X(i,l)*X(j,l) =E= 0;

eq_restrict_arrival_linear(i,j,l)$(P(i,j)).. X(i,l) + X(j,l) =L= 1;

Model fga /all - eq_restrict_arrival_linear/;

Model fga_l /all - eq_restrict_arrival/;

*option miqcp=cplex;
*Solve fga minimizing TOTCOST using miqcp;

option limrow=0, limcol=0;

$set penalty 650 !!This comes from the Scalar - pen_one_gate

$batInclude '..\qubo_solve.gms' fga_l miqcp min TOTCOST 650 -solver=cplex -timeLimit=60 -numThreads=4

display X.l, TOTCOST.l;