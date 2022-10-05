function [mean_pkts_buff,pkts_txted]=find_moments_mean_alt_cc(w,m,mbar,p,r)
global timeq;global lambda;global na_max;
n_pkts_drop=0;
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
    if(s<=m)
        if(p~=0)
            outer_term=sum*(p^s)*(1-p);%/(1-(p^(m+1)));
            sum1=sum1+outer_term;
        else
            sum1=0;
        end
    end
    del_iter=timeq*(sum);%/(1-(p^(m+1)));
    prob=((p^s)*(1-p))/(1-(p^(m+1)));
    mean_pkts_buff=mean_pkts_buff+  (prob*lambda*del_iter*(10^(-6)));
    pkts_txted=pkts_txted + (prob*find_n_txed_min(del_iter*(10^(-6)),lambda,na_max));
end
pkts_txted=min(mean_pkts_buff,na_max);


% na_fresh=max(n_pkts_suc,0);
% na_old_left=na_fresh-min(na_max,na_fresh);
% na_tot=na_fresh+na_old_left;
% if(na_tot<na_max)
%     na_theory=na_tot;
% else
%     na_theory=na_max;
% end

% 
% %%%%%%%%%%%%%%% drop delay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum1=0;
% for s=0:1:m
%     sum=0;
%     for i=0:1:s
%         if(i<mbar)
%             w_eff=(2^(i))*w;
%         else
%             w_eff=(2^(mbar))*w;
%         end
%         l=w_eff;
%         j=ceil((l-1)/r);
%         mdash=(l-1)-(r*(j-1));
%         term1=1/l;
%         term2=j*mdash/l;
%         term3=r*(j)*(j-1)/(2*l);
%         term=term1+term2+term3;
%         sum=sum+term;
%     end
%     % sum1=sum1+(sum);
%       n_pkts=find_n_pkts_accumulated2(lambda,(sum)*(10^-6)*timeq);
%     sum1=sum1+(n_pkts);%*(p^(m+1)));
% end
% % sum1=sum1*(p^m);
% n_pkts_drop=sum1*((p^(m+1)));
