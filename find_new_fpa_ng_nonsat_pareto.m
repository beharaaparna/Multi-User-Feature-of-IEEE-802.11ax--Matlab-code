function [tau,gamma]=find_new_fpa_ng_nonsat_pareto()
options = optimset('Display','off');
x0=[0.5,0.5];
fun1 = @root2d_new_fpa_ng_23_nonsat_pareto;%root2d_new_fpa1;%root2d_new_fpa_ng_2;
as= fsolve(fun1,x0,options);
% fun1 = @root2d_new_fpa_ng;
% gamma = fzero(fun1,[0.01,0.99]);
% tau2 =find_tau_fpa_ng(as(1));
tau=as(2);gamma=as(1);
end