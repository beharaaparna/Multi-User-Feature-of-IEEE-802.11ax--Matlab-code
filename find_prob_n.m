function [prob,flag_error]=find_prob_n2(n,lambda,K_scale,alpha,t)
global fn_array;
[fn,flag_error]=find_fn_integral(n,lambda,K_scale,alpha,t);
fn_array=[fn_array fn];
prob=fn/lambda;
if(prob >1)
    prob=1;
end
end