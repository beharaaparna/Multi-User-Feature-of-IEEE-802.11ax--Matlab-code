function rr=find_rr_pareto(k,alpha,time)
v1=k/(k+time);
v2=v1^alpha;
v3=1-v2;
rr=v3;
end