function indukaltmatrixnorma(A,norma)
%% Adott indukalt matrixnorma szamitas switch-case szerkezettel
%
% Hibauzenet nem negyzetes matrixokra
% Hibauzenet rosszul gepelt sztring eseten
%

if size(A,1)==size(A,2)
    switch norma
        case 'egyes' 
            norma_1=norm(A,1)
        case 'kettes'
            norma_2=norm(A,2)
        case 'vegtelen'
            norma_inf=norm(A,'inf')
        otherwise
            warning('Nem jo sztringet gepeltel. Kerlek az egyes, kettes vagy vegtelen szavak valamelyiket gepeld!')
    end
    else
    disp('Kerlek negyzetes matrixot gepelj be!')
end
