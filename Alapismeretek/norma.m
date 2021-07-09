function [norma_1, norma_2, norma_inf, norma_fro] = norma(A)

%% Bemenõ paraméterek listája

% A         matrix

%% Kimenõ paraméter

% norma_1     indukalt 1-es matrixnorma
% norma_2     indukalt 2-es matrixnorma
% norma_inf   indukalt vegtelen matrixnorma
% norma_fro   Frobenius matrixnorma

% Normak
norma_1=max(sum(abs(A),1));
norma_2=sqrt(max(eig(A'*A)));
norma_inf=max(sum(abs(A),2));
norma_fro=sqrt(sum(sum(A.^2)));



