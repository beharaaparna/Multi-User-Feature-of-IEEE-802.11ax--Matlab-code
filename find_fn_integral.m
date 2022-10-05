function [fn,flag_error]=find_fn_integral(n,lambda,beta,m,t)
delta=.1;max_val=50;flag_error=0;
sum=0;temp_array=[];
for v=0.001:delta:max_val
    val2=(t/(n*beta));
    val1=exp(-(1+val2)*v);
     if(v>2)
        here=1;
     end
    if(v==132)
        here=1;
    end
    val3=find_phi(m,n,v/n); 
     temp_array=[temp_array val3];
    if(isinf(val3))
        here=1;
    end
    val4=val3*val1;
    sum=sum+(delta*val4);
end
sumq=sum/(n*beta);
if(sumq>1)
    flag_error=1;
%     sumq=1;
end
fn=max(0,sumq);
end