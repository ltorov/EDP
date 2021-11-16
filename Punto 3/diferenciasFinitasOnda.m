
%Diferencias finitas para onda
ax = 0;
bx = 1;
at = 0;
bt = 1;
jj = 6;
c = 4;
hf = @(j) 2^(-j);
h = hf(jj);
k=(bt-at)/c;
M = (bx-ax)/h;
N = M;
uexact = @(x,t) sin(pi*x).*cos(4*pi*t);
f = @(x,t) 16*pi^2*sin(pi*x).*cos(4*pi*t);

g1=@(x) sin(pi*x);

g2= 0;

g3= 0;

g4= 0;

m=M+1;n=N+1; mn=m*n;




X=(ax:h:bx);

T=at+(0:N)*k;

A=zeros(mn,mn);
b=zeros(mn,1);

for i=2 :m-1 

for j=2:n-1

A(i+(j-1)*m,i-1+(j-1)*m)=1/h^2;
A(i+(j-1)*m,i+1+(j-1)*m)=1/h^2;

A(i+(j-1) *m,i+(j-1)*m)=-2/h^2-2/k^2;

A(i+(j-1)*m,i+(j-2)*m)=1/k^2;
A(i+(j-1)*m,i+j*m)=1/k^2;

b(i+(j-1)*m)=f(X(i),T(j));

end

end

for i=1:m % bottom and top boundary points

A(i,i)=1;
b(i)=g1(X(i));

A(i+(n-1)*m,i+(n-1)*m)=1;
b(i+(n-1)*m)=g2;

end

for j=2:n-1 % left and right boundary points

A(1+(j-1)*m,1+(j-1)*m)=1;
b(1+(j-1)*m)=g3;
A(m+(j-1)*m,m+(j-1)*m)=1;
b(m+(j-1)*m)=g4;

end

U=A\b;





Tnew = repmat(T',size(T,2),1);

Xnew = zeros(size(X,2)*size(X,2),1);

cont = 1;

for i=1:size(X,2)
    for j = 1:size(X,2)
        Xnew(cont)=X(i);
        cont = cont+1;
    end
end

Uexact = uexact(Xnew,Tnew);
Error = U-Uexact;
results = [Xnew Tnew U Uexact Error];
variablenames = {'x','t', 'U aproximada', 'U exacta', 'Error absoluto'};
results = array2table(results, 'VariableNames',variablenames);
disp(results)

tablalatex.data = results;
tablalatex.tableColLabels = variablenames;
latex = latexTable(tablalatex);

w=reshape(U(1:mn),m,n);
meshu = mesh(X,T,w)
