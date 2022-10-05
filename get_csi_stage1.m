function [time_csi] = get_csi_stage1(n,r_csi,ysc,an_ap,nsa)
global Nru;r_ul=r_csi;global K;
global time_csi;global t_breport; global t_ndpa; global t_ndp; global t_trigger_brpoll;
% n=n+(Nsa/2);
% n_stas_per_stage=(floor(n/(74)))
t_phy_leg=20;%microsec
t_sifs=16;
aifs_csi=25;
l_sf=16;%bits
l_ndpa=168+(32*n);
% l_ndpa=168;
l_tb=18;
r_leg=6;%in mbps
t_ndp=168;
l_trigger_brpoll=224+(48*n);
% l_trigger_brpoll=224;
t_ndpa=t_phy_leg+ceil((l_sf+l_ndpa+l_tb)/r_leg);
t_trigger_brpoll=t_phy_leg+ceil(l_trigger_brpoll/r_leg);
t_phy_he=228;%microsec
l_mh=320;
n_ang=56;b1=8;b2=8;Nsg=16;
l_breport=64+n_ang*((b1+b2)/2)*(ysc/Nsg);
t_breport=t_phy_he+ceil((l_sf+l_mh+l_breport+l_tb)/r_ul);
time_csi_unsolicit=t_ndpa+t_ndp+t_trigger_brpoll+t_breport;
if(n==1)
    time_csi_unsolicit=t_ndpa+t_ndp+t_breport;
end
% (floor(n/(Nru)))
K=ceil((n+nsa)/(74));
time_csi=aifs_csi + t_ndpa+t_sifs+t_ndp+t_sifs+(ceil((n+nsa)/(74)))*(t_sifs+t_trigger_brpoll+t_sifs+t_breport);
end