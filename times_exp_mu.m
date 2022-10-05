function [time_tx_mu_ul] = times_exp_mu(vu)
global r_ul;
global PayLoad;
global an_ap;
global Vs;
global an_sta;
global nra;
global max_na;
global Nsa;
global mcs;
global cur_na;
global mh;global md;
global t_tfr;global t_bsr;global t_trigger;global t_mu_data_ul;global t_msback;

payload = PayLoad;vs=Vs;
t_sifs = 16 ;% in microsec;
t_leg_pre = 20 ;%in microsec;
t_mu_he_pre=168;
sf = 16;%in bits
tb=18;% in bits
md=32;%in bits
mh=320;%in bits;
r_leg = 6;

ba_tri=224+(48*vu);
msback = 176+(288*vu);
Ltf = (50 + (10*Nsa))*8;%in bits
Lbsr=32*8;%in bits
Ltf=140*8;
t_trigger =  ceil(((ba_tri+sf+tb))/ r_ul) + t_leg_pre;
t_msback = ceil(((msback+sf+tb))/ r_ul) + t_leg_pre;
t_tfr =  ceil(((Ltf+sf+tb))/ r_leg) + t_leg_pre;
t_bsr =  ceil(((Lbsr+sf+tb))/ r_ul) + t_leg_pre;

t_trigger= 100;%micro sec
t_msback = 40;
t_tfr=70;%doubtful

[na] = find_na_mu();
cur_na=na;
if (na ==1)
    t_mu_data = ceil(((sf+na*(mh+payload)+tb))/ r_ul)  + t_mu_he_pre;
else
    t_mu_data = ceil(((sf+na*(mh+md+payload)+tb))/ r_ul)  + t_mu_he_pre;
end
t_mu_data=5300;
t_mu_data_ul=t_mu_data;
time_tx_mu_ul = t_tfr + t_sifs + t_bsr  +t_sifs +t_trigger +t_sifs +t_mu_data +t_sifs + t_msback;%in microsecs
% time_tx_mu_ul=5981;
end