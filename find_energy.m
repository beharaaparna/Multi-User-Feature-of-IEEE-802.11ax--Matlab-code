function [energy_eff,energy] = find_energy(total_stas,twt_flag,n,th,timeq,factor,t_g,n_sa,n_suc)
global total_nodes;
global p_r;
global p_t;
global act_stas;
global p_i;
global p_d;
global total_time;global ra1;
global ra;
global r_ul;
global PayLoad;
global an_ap;
global Vs;
global an_sta;
global max_na;
global Nsa;
global mcs;
global cur_na;
global mh;global md;
global t_tfr1;global t_bsr1;global t_trigger1;global t_mu_data_ul1;global t_msback1;global timeq1;
global s_throu;

if(twt_flag)
   act_stas= total_stas;
end

total_nodes=n;
payload = PayLoad;vs=Vs;
t_sifs = 16 ;% in microsec;
t_leg_pre = 20 ;%in microsec;
t_mu_he_pre=168;
sf = 16;%in bits
tb=18;% in bits
md=32;%in bits
mh=320;%in bits;
r_leg = 6;
p_s=ra;
p_r =600;p_t =1000;p_i=300;p_d=150;
% n_suc=ra;
n_tot=total_nodes;
t_mu_data=t_mu_data_ul1;
if(n~=0)
e1=n_tot*((t_tfr1*p_r)+ (4*(t_sifs*p_i)))
e1_1=t_g*((t_bsr1*p_t) + (t_trigger1*p_r))
e2= ((n_tot-n_suc)*(t_mu_data+t_msback1)*p_i)
e2_1= ((n_tot-t_g)*(t_bsr1+t_trigger1)*p_i)
e3=(n_suc*((t_mu_data*p_t)+(t_msback1*p_r)))
e4=((act_stas-n_tot-n_sa)*timeq*p_d);
e= (e1+e1_1+e2+e2_1+e3+ e4)/(act_stas);
else
    e4=((act_stas-n_tot-n_sa)*timeq*p_d)
    e=e4;
end
e1_s=(t_trigger1+t_msback1)*p_r;
e2_s=((4*t_sifs) + t_bsr1)*p_i;
e3_s=p_t*t_mu_data;
e_sa=((n_sa)/(act_stas))*(e1_s +e2_s +e3_s);
e=e+e_sa;
energy_eff =th/(n_tot*e);
energy=(e/(timeq));
