function qq=find_qq_pareto(k,alpha,time)
lam=(alpha-1)/k;
vv=lam*time;
qq=min(1,vv);
end