function [K,alpha]= find_pareto_lambda(lambda)
alpha=2;%1.01;%rand(2,10);
K=(alpha-1)/lambda;
end