function [th,ra]= func_new_mc(n)
global W_ng;
global W_g;
global Nra_g;
global m_ul;
global an_ap;
global g;
global m_limit;
global Nra_ng;
global timeq;
global Nra;
global Nsa;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r_ul;
global PayLoad;
global an_ap;
global Vs;
global an_sta;
global max_na;
global Nsa;
global mcs;
global cur_na;
global tmax;
global b;

global gamma_g_mc;
global gamma_ng_mc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_station=n;packet=PayLoad;tau_g=0;nra=n;
ap_ant=an_ap;
if(an_ap~=1)
    if(Nra_g~=0)
[tau_g,gamma_g]=find_new_mc_g();gamma_g_mc=gamma_g;
    else
        tau_g=0;gamma_g=0;
    end
end
if(Nra_ng~=0)
[tau_ng,gamma_ng]=find_new_mc_ng();gamma_ng_mc=gamma_ng;
else
    tau_ng=0;gamma_ng=0;
end
sum11=0;ag=(tau_g/Nra_g);
if(Nra_ng~=0)
    a_ng=tau_ng/Nra_ng;
    if(an_ap==1)
        a_ng=tau_ng/Nra;
    end
else
    a_ng=0;
end
if(ap_ant~=1)
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
end
sum1=0;
if(ap_ant~=1)
    if(n~=g)
        if(Nra_ng~=0)
            ps_ng_= (nchoosek(n-g,1))*(a_ng)*((1-a_ng)^(n-g-1));
        else
            ps_ng_=0;
        end
        Ps1=(Nra_ng*ps_ng_)+(Nra_g*sum11);
    else
        ps_ng=0;
        Ps1 =sum11 + ps_ng;
    end
else
    ps_ng_=nra*(a_ng)*((1-a_ng)^(nra-1));
    Ps1=(Nra*ps_ng_) ;ps_total=0;
end

ps_ana=(Ps1);
[T1] = times_exp_mu(((Nsa*an_ap)+(ps_ana)));
if(Nra_g~=0)
s_throu=(((Nsa*an_ap)+(ps_ana))*(packet));th=s_throu/timeq;
else
    s_throu=(((Nsa*an_ap)+(ps_ana))*(packet));th=s_throu/timeq;

end
s_bsr=ps_ana;ra=s_bsr/nra;
end