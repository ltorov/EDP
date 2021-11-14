
%Backward Euler
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
h = hf(2);

T = (a:h:b)';
n = size(T,1);
Y=zeros(n,N);
Y(1,:)=y0;

for j=1:n-1
    back = @(w) w - h*y(T(j+1),w)' - Y(j,:);
    B = fsolve(back,Y(j,:));
    Y(j+1,:) = B';
end

Y1exact = y1exact(T);
Y2exact = y2exact(T)

Yexact = [Y1exact Y2exact];


plot(T', Y(:,1))
hold on
plot(T', Y1exact)
legend('Solución aproximada','Solución exacta')
title("Aproximación de y1(t) con Euler implícito" )


plot(T', Y(:,2))
hold on
plot(T', Y2exact)
legend('Solución aproximada','Solución exacta')
title("Aproximación de y2(t) con Euler implícito" )