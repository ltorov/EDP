
a= 0;
b = 1;
hfun = @(j) 2.^(-j);
h = hfun(1)
X = (a:h:b)
n = size(X,2)-1
epsilon = 0.1 % este valor se varía
uexact = @(x) 1 + x + (exp(x/epsilon)-1)./(exp(1/epsilon)-1);    

ujminus1 =  1 - (h/(2*epsilon))
uj = -2
ujplus1 =  1 + (h/(2*epsilon))
A = zeros(n,n);
A = diag(uj*ones(1,n)) + diag(ujplus1*ones(1,n-1),1) + diag(ujminus1*ones(1,n-1),-1)
fj = - h^2/epsilon
C = fj
C(1) = fj -  (1 + (h/(2*epsilon)))
C(end) =  fj -  3*(1 - (h/(2*epsilon)))

U = ones(n+1,1)
U(1) = 1
U(end) = 3
U([2:end-1]) = tridiag(uj,ujminus1,ujplus1,C)

Uexact = uexact(X)'
Error = U-Uexact

results = [X' U Uexact Error];
var_names = {'x', 'U aproximada', 'U exacta', 'Error absoluto'};
dat_table = array2table(results, 'VariableNames',var_names);
disp(dat_table)

plot(X', U)
hold on
plot(X', Uexact)

legend('Solución aproximada','Solución exacta')
title("Aproximación de U(x) con j=1" )