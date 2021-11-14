
%forward euler
N= 2;
y = @(t,y) [9*y(1)+24*y(2)+5*cos(t)-1/3*sin(t); -24*y(1)-51*y(2)-9*cos(t)+1/3*sin(t)];
yexact = @(t,y) [2*exp(-3*t)-exp(-39*t)+1/3*cos(t); -exp(-3*t)+2*exp(-39*t)-1/3*cos(t)];

a = 0;
b = 1;
y10 = 4/3;
y20 = 2/3;
y0 = [y10; y20];
hf = @(j) 2.^(-j);
h = hf(2);

X = (a:h:b)';
n = size(X,1);
Y=zeros(n,N);
Y(1,:)=y0;

for j=1:n-1
    Y(j+1,:)=Y(j,:)+h*y(X(j),Y(j,:))';
end

X
Y
Yexact  = yexact(X)
