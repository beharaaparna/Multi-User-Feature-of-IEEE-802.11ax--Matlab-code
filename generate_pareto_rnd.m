function val= generate_pareto_rnd(K,alpha,n)
u=rand(1,n);
u1=u .^ (1/alpha);
u2=1 ./ u1;
u3=u2-1;
u4=K .* u3;
val=u4;
end