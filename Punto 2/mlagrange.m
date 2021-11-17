function ecuacion=mlagrange(Xi,Uaprox)
n=length(Xi);
syms x;
for i=1:n
    Li=1;
    for j=1:n
        if j~=i
            Li=Li*((x-Xi(j))/(Xi(i)-Xi(j)));
        end
    end
  L(i)=Li;
end
ecuacion=0;
for i=1:n
    ecuacion=L(i)*Uaprox(i)+ecuacion;
end
ecuacion=simplify(expand(ecuacion));

end
