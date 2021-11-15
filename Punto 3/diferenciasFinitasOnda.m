
%Diferencias finitas para onda
x1 = 0;
xr = 1;
yb = 0;
yt = 1;
M = 20;
N = 20;
xr,yb,yt,M ,N
f=@(x,t) 16*pi^2*sin(pi*x)*cos(4*pi*t); % define input function data

g1=@(x) sin(pi*x); % define boundary values

g2=@(x) 0;

g3=@(t) 0;

g4=@(t) 0;

m=M+1;n=N+1; mn=m*n;

h= (xr-x1)/M;h2=h^2;k=(yt-yb)/4;k2=k^2 ;

x=x1+(0:M)*h; % set mesh values

y=yb+(0:N)*k;

A=zeros(mn,mn);b=zeros(mn,1);

for i=2 :m-1 % interior points

for j=2:n-1

A(i+(j-1)*m,i-1+(j-1)*m)=1/h2;
A(i+(j-1)*m,i+1+(j-1)*m)=1/h2;

A(i+(j-1) *m,i+(j-1)*m)=-2/h2-2/k2;

A(i+(j-1)*m,i+(j-2)*m)=1/k2;
A(i+(j-1)*m,i+j*m)=1/k2;

b(i+(j-1)*m)=f(x(i),y(j));

end

end

for i=1:m % bottom and top boundary points

j=1; A(i+(j-1)*m,i+(j-1)*m)=1;b(i+(j-1)*m)=g1(x(i));

j=n;A(i+(j-1)*m,i+(j-1)*m)=1;b(i+(j-1)*m)=g2(x(i));

end

for j=2:n-1 % left and right boundary points

i=1;A(i+(j-1)*m,i+(j-1)*m)=1;b(i+(j-1)*m)=g3(y(j));

i=m;A(i+(j-1)*m,i+(j-1)*m)=1;b(i+(j-1)*m)=g4(y(j));

end

v=A\b;% solve for solution in v labeling

w=reshape(v(1:mn),m,n); %translate from v to w

mesh(x,y,w)

