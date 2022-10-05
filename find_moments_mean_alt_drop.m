function [mea_del,mean_delay_drop]=find_moments_mean_alt(w,m,mbar,p,r)
global timeq;
sum1=0;%find_second_moments
% m=m-1;
for s=0:1:m
    sum=0;
    for i=0:1:s
        if(i<mbar)
            w_eff=(2^(i))*w;
        else
            w_eff=(2^(mbar))*w;
        end
        l=w_eff;
        j=ceil((l-1)/r);
        mdash=(l-1)-(r*(j-1));
        term1=1/l;
        term2=j*mdash/l;
        term3=r*(j)*(j-1)/(2*l);
        term=term1+term2+term3;
        sum=sum+term;
    end
    if(s<=(m))
        if(p~=0)
            outer_term=sum*(p^s)*(1-p);
            sum1=sum1+outer_term;
        else
            sum1=0;
        end
    else
        gh=(1-(p^(m+1)));
%         outer_term=sum*(p^s)*(1-p)/gh;
%         sum1=sum1+outer_term;
%         wdd=p^(m+1);
%         wdash=1-wdd;
%         outer_term=sum*(p^s)*(1-p)/(1-(p^(m+1)));
%         sum1=sum1+outer_term;
    end
end
mea_del=timeq*sum1/(1-(p^(m+1)));

%%%%%%%%%%%%%%% drop delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum1=0;
for s=0:1:m
sum=0;
for i=0:1:s
    if(i<mbar)
        w_eff=(2^(i))*w;
    else
        w_eff=(2^(mbar))*w;
    end
    l=w_eff;
    j=ceil((l-1)/r);
    mdash=(l-1)-(r*(j-1));
    term1=1/l;
    term2=j*mdash/l;
    term3=r*(j)*(j-1)/(2*l);
    term=term1+term2+term3;
    sum=sum+term;
end
sum1=sum1+(sum);
% sum1=sum;%+(sum*(p^s));
end
% sum1=sum1*(p^m);
mean_delay_drop=timeq*sum1*((p^(m+1)));
