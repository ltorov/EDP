%Forward Euler
N= 2;
y = @(t,y) [9*y(1)+24*y(2)+5*cos(t)-1/3*sin(t); -24*y(1)-51*y(2)-9*cos(t)+1/3*sin(t)];
yexact = @(t) [2*exp(-3*t)-exp(-39*t)+1/3*cos(t); -exp(-3*t)+2*exp(-39*t)-1/3*cos(t)];
y1exact = @(t) 2*exp(-3*t)-exp(-39*t)+1/3*cos(t);
y2exact = @(t) -exp(-3*t)+2*exp(-39*t)-1/3*cos(t);

a = 0;
b = 1;
y10 = 4/3;
y20 = 2/3;
y0 = [y10; y20];
hf = @(j) 2.^(-j);
for jj=6:8
    h = hf(jj);
    
    T = (a:h:b)';
    n = size(T,1);
    Y=zeros(n,N);
    Y(1,:)=y0;
    
    for j=1:n-1
        Y(j+1,:)=Y(j,:)+h*y(T(j),Y(j,:))';
    end
    plot(T', Y(:,2))
    hold on
    
    Y1exact = y1exact(T);
    Y2exact = y2exact(T);

    Error1 = abs(Y1exact-Y(:,1));
    Error2 = abs(Y2exact-Y(:,2));
 

    results = [T Y(:,1) Y(:,2) Y1exact Y2exact Error1 Error2];
    variablenames = {'t', 'y1 aproximada','y2 aproximada','y1 exacta','y2 exacta', 'Error absoluto y1','Error absoluto y2'};
    results = array2table(results, 'VariableNames',variablenames);
    disp(results);

    tablalatex.data = results;
    tablalatex.tableColLabels = variablenames;
    latex = latexTable(tablalatex);
    
end

Y1exact = y1exact(T);
Y2exact = y2exact(T);

plot(T', Y2exact)
legend('Solución aproximada con h=0.0312',...
    'Solución aproximada con h=0.0156',...
    'Solución aproximada con h=0.0078', ...
    'Solución aproximada con h=0.0039', ...
    'Solución aproximada con h=0.0020', ...
    'Solución exacta')
title("Aproximación de y2(t) con Euler explícito" )





