
% Adams-Bashforth de 3 pasos
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
for jj=6:9
    h = hf(jj);
    
    T = (a:h:b)';
    n = size(T,1);
    Y=zeros(n,N);
    Y(1,:) = y0;
    
    %Uso Runge-Kutta para encontrar los primeros términos
    for j=1:3
            k1 = y(T(j),Y(j,:))';
            k2 = y(T(j)+C(2)*h,Y(j,:)+A(2,1)*k1*h)';
            k3 = y(T(j)+C(3)*h,Y(j,:)+(A(3,1)*k1+A(3,2)*k2)*h)';
            k4 = y(T(j)+C(4)*h,Y(j,:)+(A(4,1)*k1+A(4,2)*k2+A(4,3)*k3)*h)';
            k5 = y(T(j)+C(5)*h,Y(j,:)+(A(5,1)*k1+A(5,2)*k2+A(5,3)*k3+A(5,4)*k4)*h)';
            Y(j+1,:)=Y(j,:)+h*(B(1)*k1+B(2)*k2+B(3)*k3+B(4)*k4+B(5)*k5);
    end
    
    
    for j=3:n-1
        Y(j+1,:)=Y(j,:)+h*((5/12)*y(T(j-2),Y(j-2,:))' - (4/3)*y(T(j-1),Y(j-1,:))' + (23/12)*y(T(j),Y(j,:))');
    end
    plot(T', Y(:,2))
    hold on

end

Y1exact = y1exact(T);
Y2exact = y2exact(T);

plot(T', Y2exact)
legend('Solución aproximada con h=0.0156',...
    'Solución aproximada con h=0.0078',...
    'Solución aproximada con h=0.0078',...
    'Solución aproximada con h=0.0020',...
    'Solución exacta')
title("Aproximación de y2(t) con Adams Bashforth de 3 pasos" )

