function x=newtonrendszer(f,x0,maxit,TOL,fJ);

%% Input
% fJ  - a fuggveny pontos Jacobija (opcio, egyebkent Numerikus Jacobi)
    
numJ=nargin<5;
if nargin<3, maxit=100; end;
if nargin<4, TOL=1e-6; end;

x=x0; i=0; dx=ones(size(x0));   
while norm(f(x))>TOL && norm(dx)>TOL && i<maxit
  if numJ, J=NumerikusJacobi(f,x); else J=fJ(x); end;
  dx=-J\f(x); x=x+dx; i=i+1;
end
if i>=maxit 
  error('A Newton modszer nem konvergal: a maximalis iteracioszamot elertuk'); 
end;
