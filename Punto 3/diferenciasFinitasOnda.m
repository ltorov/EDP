
%Diferencias finitas para onda
ax = 0;
bx = 1;
at = 0;
bt = 1;
M = 20;
N = 20;
f=@(x,t) 16*pi^2*sin(pi*x)*cos(4*pi*t); % define input function data

g1=@(x) sin(pi*x); % define boundary values

g2= 0;

g3= 0;

g4= 0;

m=M+1;n=N+1; mn=m*n;


h= (bx-ax)/M;
k=(bt-at)/4;

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

v=A\b;% solve for solution in v labeling

w=reshape(v(1:mn),m,n); %translate from v to w

mesh(X,T,w)

