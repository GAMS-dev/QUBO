$if not set prod_of_primes $set prod_of_primes 3071
$if not set bit_length_1   $set bit_length_1      6
$eval bit_length_2 (floor(log2(%prod_of_primes%))+1)-%bit_length_1%+1

Set i           'indices for binary encoding factor 1' /0*%bit_length_1%/,
    j           'indices for binary encoding factor 2' /0*%bit_length_2%/;

Parameters
    coeff_x(i) 'coefficients of factor 1',
    coeff_y(j) 'coefficients of factor 2';
    
coeff_x(i) = power(2,ord(i)-1);
coeff_y(j) = power(2,ord(j)-1);

Binary Variable X(i), Y(j);

Free Variable TOTCOST;

Equation
    cost_fn 'Objective function',
    obj;

cost_fn.. TOTCOST =E= 0;
obj.. sum(i, coeff_x(i)*X(i)) * sum(j, coeff_y(j)*Y(j)) =E= %prod_of_primes% ;

Model pf /all/;

option MIQCP=gurobi;
Solve pf using MIQCP min TOTCOST;

Parameter report(*);

report('Factor 1') = sum(i,coeff_x(i)*X.l(i));
report('Factor 2') = sum(j,coeff_y(j)*Y.l(j));

option report:2:0:1; display report;