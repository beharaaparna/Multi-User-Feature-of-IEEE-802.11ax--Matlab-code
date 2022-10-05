function [mea_del,pkts_txted]=find_moments_mean_alt_edited1(w,m,mbar,p,r)
global timeq;global lambda;global na_max;global K_scale;global alpha;
sum1=0;%find_second_moments
% m=m-1;
mean_pkts_buff=0;pkts_txted=0;
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
    end
    del_iter=timeq*(sum);%/(1-(p^(m+1)));
    prob=((p^s)*(1-p))/(1-(p^(m+1)));
%     mean_pkts_buff=mean_pkts_buff+ (prob*find_avg_pkts_pareto(lambda,K_scale,alpha,del_iter*(10^(-6))));%lambda*del_iter*(10^(-6)));
    pkts_txted=pkts_txted + (prob*find_n_txed_min(del_iter*(10^(-6)),lambda,na_max))
end
% if(pkts_txted>mean_pkts_buff)
%     pkts_txted=mean_pkts_buff;
% end
if(pkts_txted>na_max)
    pkts_txted=na_max;
end
mea_del=timeq*sum1/(1-(p^(m+1)));
end


% %%%%%%%%%%%%%%% drop delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum1=0;
% for s=0:1:m
% sum=0;
% for i=0:1:s
%     if(i<mbar)
%         w_eff=(2^(i))*w;
%     else
%         w_eff=(2^(mbar))*w;
%     end
%     l=w_eff;
%     j=ceil((l-1)/r);
%     mdash=(l-1)-(r*(j-1));
%     term1=1/l;
%     term2=j*mdash/l;
%     term3=r*(j)*(j-1)/(2*l);
%     term=term1+term2+term3;
%     sum=sum+term;
% end
% sum1=sum1+(sum);
% % sum1=sum;%+(sum*(p^s));
% end
% % sum1=sum1*(p^m);
% mean_delay_drop=timeq*sum1*((p^(m+1)));
