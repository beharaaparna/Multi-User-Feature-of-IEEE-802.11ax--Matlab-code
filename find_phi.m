function phi=find_phi(m,n,v)
val1=(-1)^(n+1);pi=3.14;
val2=m^n;
lim=floor((n-1)/2);
sum=0;
for r=0:1:lim
    a1=(-pi^2)^r;
    a2=nchoosek(n, (2*r)+1);
    a3=find_ei(m+1,v);
    if(isinf(a3))
        break;
    end
    a4=a3^(n-(2*r)-1);
    a5_num=v^m;a5_den=factorial(m);
    a5=a5_num/a5_den;
    a6=a5^((2*r)+1);
    full=a1*a2*a4*a6;
    sum=sum+full;
end
phi=val1*val2*sum;
end