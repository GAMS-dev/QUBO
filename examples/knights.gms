$title Maximum Knights Problem (KNIGHTS,SEQ=158)

$onText
This MIP model finds the maximum number of knights that can be
placed on a board. Two different formulations are presented.
The second formulation is 'tight' and may perform better with certain
MIP codes. Once we found the max number of knights, we solve a series
of MIPs to find ALL solutions.

We will use lags (relative positions) to describe the allowed moves.
The labels H and V indicate horizontal and vertical moves as shown
below:

                 0 0
                0   0
                  X
                0   0
                 0 0


Dudeney, H E, Amusements in Mathematics. Dover, New York, 1970.

Keywords: mixed integer linear programming, maximum knights problem, mathematics
$offText

Set
   i 'size of board'            /  1*8  /
   n 'number of possible moves' / m1*m8 /;

Alias (i,j,k);

Table move(*,n) 'all possible knight moves'
      m1 m2 m3 m4 m5 m6 m7 m8
   H  -2 -2 -1 -1 +1 +1 +2 +2
   V  -1 +1 -2 +2 -2 +2 -1 +1;

Variable total;

Binary Variable x(i,j);

Equation
   deftotal        'total knights on board'
   defmove(i,j)    'move restrictions'
   defmovex(n,i,j) 'move restrictions';

deftotal..        total =e= sum((i,j), x(i,j));

defmove(i,j)..    sum(n, x(i + move('h',n),j + move('v',n))) =l= card(i)*(1 - x(i,j));

defmovex(n,i,j).. x(i + move('h',n),j + move('v',n)) =l= 1 - x(i,j);

Model
   knight  / deftotal, defmove  /
   knightx / deftotal, defmovex /;

option optCr = 0, optCa = .999;

*solve knight using mip max total;

option limrow=0, limcol=0;

$batinclude  '..\qubo_solve.gms' knight mip max total 1 -method=qpu -timeLimit=10

display total.l;