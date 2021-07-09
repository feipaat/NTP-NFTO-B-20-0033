function [k x]=newton(fuggveny,x,maxit,TOL,derivalt)

format long
hiba=TOL+1;
k=0;
while k<maxit && hiba>TOL
   k=k+1;
   xregi=x;
   x=x-eval(fuggveny)/eval(derivalt);
   hiba=abs(x-xregi);
end;

