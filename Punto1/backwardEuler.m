a = 0;
b = 2;
f = @(x,z) z./(x+1)-x.^2+1;
y0 = 1;
hf = @(j) 2.^(-j);
h = hf(1);
X = a:h:b;
M = size(X,2);
Y = zeros(M,1);
Y(1) = y0;

for j=1:M-1
    back = @(w) w - h*f(X(j+1),w) - Y(j);
    xj=X(j+1);
    Y(j+1)=fsolve(back,Y(j));
end

[X' Y]







