function na_theory_min=find_n_txed_min(delay,lambda,na_max)
sum1=0;global K_scale;global alpha;sum_prob=0;prob_array=[];
for i=0:1:na_max-1
    [valo1,err_flag]=find_prob_n(i,lambda,K_scale,alpha,delay);
    if(valo1<0)
        break;
    end

    if(isnan(valo1))
        break;
    end

      if(valo1>1)
        break;
    end
    prob_n=valo1;%find_pois_prob_n_pkts(i,lambda,delay);
    sum_prob=sum_prob+prob_n;
    prob_array=[prob_array valo1];
    if(sum_prob>1)
        break;
    end
    if(isnan(prob_n))
        break;
    end
         if(isinf(prob_n) )
            break; 
    end
    sum1=sum1+(i*prob_n);
end

% sum2_prob=0;
% for j=0:1:na_max-1
%       [vv,err_flag]=find_prob_n(j,lambda,K_scale,alpha,delay);
% %     vv=find_pois_prob_n_pkts(j,lambda,delay);
%     if(isnan(vv) )
%             break; 
%     end
%      if(isinf(vv) )
%             break; 
%     end
%     sum2_prob = sum2_prob  + vv;
%     
% end
if(sum_prob>1)
    sum_prob=1;
end

sum2=(na_max*(1-sum_prob));

na_theory_min=(sum1)+sum2;
if(isnan(na_theory_min))
    here=1;
end
end