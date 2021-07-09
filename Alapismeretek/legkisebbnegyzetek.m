function [x,Jnorma,MatlabPolyfit]=legkisebbnegyzetek(t,f,n)
%% A Gauss-féle normál-egyenlet megvalósítása n-edfokú polinomok esetén

% A feladat: adott (t_i,f_i) alakú pontokra legkisebb négyzetek módszerének
% megvalósítása Gauss-féle normál-egyenlettel
%
%


%% Bemenõ paraméterek listája: 

% t           alappontok
% f           t-hez tartozó értékek
% n           a közelítõ módszer rendje

%% Elõkészületek 

e=ones(max(size(t)),1);
A(:,1)=e;
A(:,2)=t';
for i=3:n+1
    A(:,i)=A(:,i-1).*t';
end
ATranA=(A'*A);
ATranf=(A'*f');
x=(A'*A)\(A'*f');
J=A*x-f';
Jnorma=norm(J,2)^2;

%% Biztonsági összevetés a Matlab beépített függvényével

MatlabPolyfit=rot90(polyfit(t,f,n)',2);


