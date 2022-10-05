function prob_n=find_pois_prob_n_pkts(i,lambda,delay)
t1=exp(-lambda*delay);
t2=(lambda*delay)^i;
t3=factorial(i);
prob_n=t1*t2/t3;
end