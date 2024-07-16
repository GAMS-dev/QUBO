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

$batinclude qubo_solve knight mip max total 1 qpu "" 10

display total.l;
*solve knightx using mip max total;

$exit

* Now we try to see how many different ways are there to arrange
* the same number of knights.

Scalar maxknight;
maxknight = total.l;
total.lo  = total.l - 0.5;

option optCr = 1, optCa = 100, limCol = 0, limRow = 0, solPrint = off;

Set
   ss    'max number of solutions groups' / seq1*seq20 /
   s(ss) 'dynamic set for solution groups';

Parameter cutval 'all possible solutions for cut generation';

Equation cut(ss) 'known solutions to be eliminated';

cut(s).. sum((i,j), cutval(s,i,j)*x(i,j)) =l= maxknight - 1;

Model knight1 / deftotal, defmovex, cut /;

s(ss) = no;
total.lo = maxknight - .5;
knight1.solveStat = %solveStat.normalCompletion%;
knight1.modelStat = %modelStat.optimal%;

loop(ss$(knight1.solveStat = %solveStat.normalCompletion% and
        (knight1.modelStat = %modelStat.optimal% or
         knight1.modelStat = %modelStat.integerSolution%)),
   s(ss) = yes;
   cutval(ss,i,j) = x.l(i,j);
   solve knight1 maximizing total using mip;
);

option  cutval:0:0:1;
display cutval;
