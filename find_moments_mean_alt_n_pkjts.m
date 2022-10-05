function [mea_del,mean_pkts_buff,pkts_txted]=find_moments_mean_alt_n_pkjts(w,m,mbar,p,r,avg_del)
global timeq;global lambda;global na_max;
mean_pkts_buff=0;pkts_txted=0;mea_del=0;
sum1=0;%find_second_moments
prob13=0;
for s=1:1:m
    prob13=prob13+((p^s)*(1-p))/(1-(p^(m+1)));
end
prob13=(1-(p^(m+1)));
old_gene=prob13*lambda*avg_del;
old_left=old_gene-min(old_gene,na_max);
if(old_left>0)
 here=1;
end

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
    end
    del_iter=timeq*(sum);%/(1-(p^(m+1)));
    prob=((p^s)*(1-p))/(1-(p^(m+1)));
    if(prob>1)
        here=1;
    end
    if((1-(p^(m+1)))>1)
        here=1;
    end

    mean_pkts_buff1=lambda*del_iter*(10^(-6)); 
    mean_pkts_buff=mean_pkts_buff + (prob*mean_pkts_buff1);
    nnn=min(mean_pkts_buff1,na_max);
%   
%        [n_tx,n_lef]=find_n_txed_min(del_iter*(10^(-6)),lambda,na_max,0);
     pkts_txted=pkts_txted+ (prob*nnn);
    
end
% pkts_txted=(prob13*min(mean_pkts_buff,na_max))+((1-prob13)*mean_pkts_buff);
old_gene_kk=mean_pkts_buff-min(mean_pkts_buff,na_max);
% pkts_txted=min(mean_pkts_buff,na_max);

% n_old=mean_pkts_buff-pkts_txted;
% n_buf=n_old+mean_pkts_buff;
% n_tx=min(n_buf,na_max);
% if(pkts_txted>mean_pkts_buff)
%     pkts_txted=mean_pkts_buff;
% end
% if(pkts_txted>na_max)
%     pkts_txted=na_max;
% end
% mea_del=timeq*sum1/(1-(p^(m+1)));