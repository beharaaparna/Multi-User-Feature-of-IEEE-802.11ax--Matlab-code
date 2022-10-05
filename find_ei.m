function term=find_ei(m,x)
if(x>65)
    here=1;
end
gamma=0.577;max_val=50;
val1_num=x^(m-1);val1_den=factorial(m-1);
val1=val1_num/val1_den;
sum=0;
for r=1:1:(m-1)
    vv=1/r;
    sum=sum+(1/r);
end
inner=gamma+log(x)-sum;
first=inner*val1;

sum1=0;
for r=0:1:max_val
    if(r~=(m-1))
        num=x^r;
        den=(r-m+1)*factorial(r);
        if(isinf(den))
            sum1=sum1;break;
        else
        ff=num/den;
        sum1=sum1+ff;
        end
    end
end
sec=sum1;
term=first+sec;
end