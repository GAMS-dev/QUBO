set facility    /chicago, boston, denver/;
    
alias (facility, i, j, k, l);

parameter
    flow(i,j)
    uflow(i,j) /chicago.boston 5,
               chicago.denver 2,
               boston.denver 3/
               
    dist(k,l)
    udist(k,l) /chicago.boston 8,
               chicago.denver 15,
               boston.denver 13/;
             
flow(i,j) = max(uflow(i,j),uflow(j,i));
dist(k,l) = max(udist(k,l),udist(l,k));

binary variable x(i,j);

free variable z;

equations obj, c1, c2;

obj.. sum((i,j,k,l), flow(i,j)*x(i,k)*x(j,l)*dist(k,l)) =E= z; 

c1(j).. sum(i, x(i,j)) =E= 1;

c2(i).. sum(j, x(i,j)) =E= 1; 

Model qap /all/;

option limrow=0, limcol=0;

$batinclude '..\qubo_solve.gms' qap MIQCP min z 200

display x.l, z.l;