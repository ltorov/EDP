a= 0;
b = 1;
epsilons = [0.1 0.01 0.001 0.0001];

epsilon = 0.0001

for jj=1:6
    hfun = @(j) 2.^(-j);
    h = hfun(jj);
    X = (a:h:b);
    n = size(X,2)-1;
    uexact = @(x) x; %cambiar esto  
    sigma = %definir que es sigma
    ujminus1 =  sigma^2;
    uj = 2-2*sigma^2;
    ujplus1 =  sigma^2;
    fj = ; %definir que es fj
    C = fj;
    C(1) = fj -  ; %definir
    C(end) =  fj -  ; %definir
    
    U = ones(n+1,1);
    U(1) = ;%definir
    U(end) = ;%definir
    U([2:end-1]) = tridiag(uj,ujminus1,ujplus1,C);
    
    Uexact = uexact(X)';
    Error = U-Uexact;
    
    results = [X' U Uexact Error];
    variablenames = {'x', 'U aproximada', 'U exacta', 'Error absoluto'};
    results = array2table(results, 'VariableNames',variablenames);
    disp(results)

    tablalatex.data = [X' U Uexact Error];
    tablalatex.tableColLabels = variablenames
    latex = latexTable(tablalatex);
    
    plot(X', U)
    hold on
    plot(X', Uexact)
    
    legend('Solución aproximada con h=0.5', ...
    'Solución aproximada con h=0.25', ...
    'Solución aproximada con h=0.125', ...
    'Solución aproximada con h=0.0625',...
    'Solución aproximada con h=0.0312',...
    'Solución aproximada con h=0.0156',...
    'Solución exacta')
    title("Aproximación de U(x) con epsilon =0.0001" )
end
