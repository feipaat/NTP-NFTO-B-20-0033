function kvadratura(a,b,m,fv,modszer)

%% A tanult harom numerikus integralo osszetett formula: erinto, trapez es Simspon
%% Input
%  a          intervallum kezdete
%  b          intervallum vege
%  m          reszintervallumok szama (Simspon m/2-re lesz alkalmazva)
%  fv         az integralando fuggveny
%  modszer    a kivant modszer


format long
h=(b-a)/m;

switch modszer
    case{'erinto'}
    	fprintf('\n');
    	X=['A feladat megoldasa ', num2str(m), ' reszintervallumon osszetett erintoformulaval tortent.'];
    	disp(X)

    	x = [a:h/2:b];
    	y = eval(fv);
    	h*sum(y(2:2:2*m))
    
    case{'trapez'}
    	fprintf('\n');
    	X=['A feladat megoldasa ', num2str(m), ' reszintervallumon osszetett trapezformulaval tortent.'];
    	disp(X)

    	x = [a:h:b];
    	y = eval(fv);
        (h/2)*(y(1)+2*sum(y(2:1:m))+y(m+1))


    
    case{'Simpson'}
	fprintf('\n');
    	X=['A feladat megoldasa ', num2str(m/2), ' reszintervallumon osszetett Simpson-formulaval tortent.'];
    	disp(X)

    	x = [a:h:b];
    	y = eval(fv);
    	(h/3)*(y(1)+4*sum(y(2:2:m))+2*sum(y(3:2:m-1))+y(m+1))

    
otherwise
        fprintf('\n');
        disp('Nem megfelelo modszer nevet gepelted be.')
end

