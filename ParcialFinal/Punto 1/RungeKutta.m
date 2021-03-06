%Runge-Kutta-Fehlberg 45 de quinto orden
N= 2;
y = @(t,y) [9*y(1)+24*y(2)+5*cos(t)-1/3*sin(t); -24*y(1)-51*y(2)-9*cos(t)+1/3*sin(t)];
yexact = @(t) [2*exp(-3*t)-exp(-39*t)+1/3*cos(t); -exp(-3*t)+2*exp(-39*t)-1/3*cos(t)];
y1exact =  @(t) 2*exp(-3*t)-exp(-39*t)+1/3*cos(t)
y2exact = @(t) -exp(-3*t)+2*exp(-39*t)-1/3*cos(t);

a = 0;
b = 1;
y10 = 4/3;
y20 = 2/3;
y0 = [y10; y20];
hf = @(j) 2.^(-j);


A = [0 0 0 0 0 0; 1/4 0 0 0 0 0;3/32 9/32 0 0 0 0;
    1932/2197 -7200/2197 7296/2197 0 0 0; 439/216 -8 3680/513 -845/4104 0 0;
    -8/27 2 -3544/2565 1859/4104 -11/40 0]
B = [16/135 0 6656/12825 28561/56430 -9/50 2/55]
C = [0 1/4 3/8 12/13 1 1/2]

for jj=4:8
    h = hf(jj);
    T = (a:h:b)';
    n = size(T,1);
    Y=zeros(n,N);
    Y(1,:)=y0;
    
    for j=1:n-1
        k1 = y(T(j),Y(j,:))';
        k2 = y(T(j)+C(2)*h,Y(j,:)+A(2,1)*k1*h)';
        k3 = y(T(j)+C(3)*h,Y(j,:)+(A(3,1)*k1+A(3,2)*k2)*h)';
        k4 = y(T(j)+C(4)*h,Y(j,:)+(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3)*h)';
        k5 = y(T(j)+C(5)*h,Y(j,:)+(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4)*h)';
        Y(j+1,:)=Y(j,:)+h*(B(1)*k1+B(2)*k2+B(3)*k3+B(4)*k4+B(5)*k5);
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
legend('Solución aproximada con h=0.0625',...
    'Solución aproximada con h=0.0312',...
    'Solución aproximada con h=0.0156',...
    'Solución aproximada con h=0.0078',...
    'Solución aproximada con h=0.0039',...
    'Solución exacta')
title("Aproximación de y2(t) con Runge Kutta Fehlberg 45 de quinto orden" )
