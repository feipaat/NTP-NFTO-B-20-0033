function Lagrange(a,b,N)
%% Programoz�si beadand� a Lagrange-interpol�ci�s polinom alkalmaz�s�ra

% A feladat:

%  Tekints�k az f(x)=1/(1+25x^2), x\in[-1,1] f�ggv�nyt. Sz�m�t�g�pen 
%  ellen�rizz�k, hogy hi�ba s�r�tj�k a feloszt�st egyenletesen a [-1,1]-en,
%  a kapott polinomsorozat nem konverg�l $f$-hez!


%% Bemen� param�terek list�ja: 

% a           intervallum kezdete
% b           intervallum v�ge
% N           interpol�ci�s polinom foka (N+1=alappontok sz�ma)


%% El�k�sz�letek

     h_N=(b-a)/(N);
     X_N=a:h_N:b;
     Y_N=1./(1+25*X_N.^2);
     
     h_2N=(b-a)/(2*N);
     X_2N=a:h_2N:b;
     Y_2N=1./(1+25*X_2N.^2);
     
     h_3N=(b-a)/(3*N);
     X_3N=a:h_3N:b;
     Y_3N=1./(1+25*X_3N.^2);
     
%% Lagrange int.pol. pol. �s f(x) �sszerak�sa

     [P_N]=lagrangepoly(X_N,Y_N);
     [P_2N]=lagrangepoly(X_2N,Y_2N);
     [P_3N]=lagrangepoly(X_3N,Y_3N);
     xx=a:h_2N^4:b;
     yy=1./(1+25*xx.^2);
     
%% Plottol�s

     plot(xx,yy,'k',xx,polyval(P_N,xx),'b',xx,polyval(P_2N,xx),'r',xx,polyval(P_3N,xx),'g');
     legend('f(x)', sprintf('L_{%g}(x)', N+1), sprintf('L_{%g}(x) ', 2*N+1), sprintf('L_{%g}(x)', 3*N+1));
     %grid
     %title('Lagrange interpol�ci�s polinomok �s az eredeti f(x) f�ggv�ny')
     axis([-0.05+a b+0.05 -0.2+a b+0.2])
     figure
     plot(xx,yy,'k',X_N,Y_N,'.b', X_2N,Y_2N,'+r',X_3N,Y_3N,'sg',...               
                'MarkerSize',8);
     legend('f(x)',sprintf('L_{%g}(x) alappontjai', N+1), sprintf('L_{%g}(x) alappontjai', 2*N+1), sprintf('L_{%g}(x) alappontjai', 3*N+1));
     %grid
     %title('Lagrange interpol�ci�s polinomok alappontjai f(x)-en')
     %axis([-0.05+a b+0.05 -0.2+a b+0.2])
     axis('square','equal')
