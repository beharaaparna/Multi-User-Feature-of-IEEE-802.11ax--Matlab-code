function [th,ra,delay1,gamma_ng,pkts_txted]= func_new_fpa_nonsat_pareto(m_limit,n,n_sa_stas)
global W_ng;
global W_g;global lambda;
global s_throu;
global Nra_g;global ra;
global m_ul;
global an_ap;
global g;
% global m_limit;
global Nra_ng;
global Nra;
global Nsa;
global timeq;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r_ul;
global PayLoad;
global an_ap;
global Vs;
global an_sta;
% global nra;
global max_na;
global Nsa;
global mcs;
global cur_na;
global tmax;
global b;
global s_throu_g;global s_throu_ng;
global gamma_ga;global gamma_nga;
global tau_nga;global tau_ga;global t_g;
global new_stas;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nra=n;
new_stas=n;
num_station=n;
packet=PayLoad;
% m_limit=m_limit-1;
ap_ant=an_ap;tau_g=0;gamma_g=0;
if(Nra_ng~=0)
    [tau_ng,gamma_ng]=find_new_fpa_ng_nonsat_pareto();
else
        tau_ng=0;gamma_ng=0;

end
% [delay1,delay2]=find_moments_mean_alt_drop(W_ng,m_limit,m_ul,gamma_ng,Nra_ng);
[delay1,pkts_txted]=find_moments_mean_alt_edited1(W_ng,m_limit,m_ul,gamma_ng,Nra_ng);
delay2=0;
tot_delay=(delay1)/2;
tau_nga=tau_ng;
if(an_ap~=1)
    if(Nra_g~=0)
    [tau_g,gamma_g]=find_new_fpa_g();
    tau_ga=tau_g;
    else
        tau_ga=0;
    end
end
 gamma_ga= gamma_g; gamma_nga= gamma_ng;
 t_g=g*tau_g;t_ng=(num_station-g)*tau_ng;
t_g=t_g+t_ng;
sum11=0;ag=(tau_g/Nra_g);
if(Nra_ng~=0)
    a_ng=tau_ng/Nra_ng;
    if(an_ap==1)
        a_ng=tau_ng/Nra;
    end
else
    a_ng=0;
end
if(g >= an_ap)
    for i=1:1:an_ap
        ps_g = (nchoosek(g,i))*(ag^i)*((1-ag)^(g-i));
        sum11=sum11+(i*ps_g);
    end
else
    for i=1:1:g
        ps_g = (nchoosek(g,i))*(ag^i)*((1-ag)^(g-i));
        sum11=sum11+(i*ps_g);
    end
end
sum1=0;
if(ap_ant~=1)
    if(nra~=g)
        if(Nra_ng~=0)
            ps_ng_= (nchoosek(nra-g,1))*(a_ng)*((1-a_ng)^(nra-g-1));
        else
            ps_ng_=0;
        end
        Ps1=(Nra_ng*ps_ng_)+(Nra_g*sum11);
    else
        ps_ng=0;ps_ng_=0;
        if(Nra_g~=0)
             Ps1 =(Nra_g*sum11);
        else
        Ps1 =sum11 + ps_ng;
        end
    end
else
    %%%%%%%% here
    qq=1;%1-min(1,lambda*tot_delay*(10^(-6)));
    nra=nra*qq;
    ps_ng_=nra*(a_ng)*((1-a_ng)^(nra-1));
    Ps1=(Nra*ps_ng_) ;ps_total=0;
end

ps_ana=(Ps1);

s_throu_g=((Nra_g*sum11))*(packet)/timeq;
s_throu_ng=((Nra_ng*ps_ng_))*(packet)/timeq;

% [T1] = times_exp_mu((((Nra+Nsa)*an_ap)));
if(Nra_g~=0)
s_throu=(((n_sa_stas)+(ps_ana))*(packet));

else
    s_throu=(((n_sa_stas)+(ps_ana))*(packet));
end
th=s_throu/timeq;
s_bsr=ps_ana;ra=s_bsr/nra;
dd=cur_na*s_throu
end