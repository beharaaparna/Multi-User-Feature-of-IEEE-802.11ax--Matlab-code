clear all;
clc;
close all;

%% Global Variables Initialisation

global r;
global act_stas;
global gamma_ga;global gamma_nga;
global PayLoad;
global W_ul;global W_g;global W_ng;
global Nra_a;global Nsa;global b;
global nra;global Nra_ng;
global m;global an_ap;global Nru;global g;global Nra_g;global an_ap;global Vs;global an_st;global max_na;
global mcs;global cur_na;
global tmax;global m_ul;global r_ul;global Nra;global m_limit;global timeq;
global gamma_g_mc;
global gamma_ng_mc;
global tau_nga;global tau_ga;
global t_g;
global t_tfr;global t_bsr;global t_trigger;global t_mu_data_ul;global t_msback;
global t_tfr1;global t_bsr1;global t_trigger1;global t_mu_data_ul1;global t_msback1;

%% Simulation Parameters given in Table.1
csi_flag=0; %Flag to Include Channel sounding
twt_flag=1;% Flag to include Target Wake Up Time
T_max=25;% TWT t_Mean
flag_mimo=0;%Flag to include MU-MIMO
lcsi=5;% Channel sounding Interval
ap_ant=1;%Number of Antennas
Nru=16;%Number of Resource Units
Nra=ceil(0.5*Nru);%Number of RA RUs
Nsa=Nru-Nra;%Number of SA RUs
grouping_ratio=0.05;% percentage of STAs selected for MU-MIMO
lim=10000;
tf_cycles =10000;% Number of time-slots (TF-Cycles)
n_iter =20;%Number of iterations
kmax=2000;
kin1=10;
k1kk=80:100:kmax;k1=k1kk;
m_limit=6;% Maximum Retry Limit
PayLoad = 12000;%payload in bits
W=64;% Minimum Contention Window Size
na=256;% Maximum Number of packets aggregated
m=5;m_ul=m;
vs=1;% Number fo antennas at the STA
mcs=11;%Modulation & Coding Index
b=160;% Bandwidth in MHz

%% Flag for MU-MIMO
if(flag_mimo)
    Nra_g=ceil(Nra*grouping_ratio);% select RUs for MU-MIMO operation
    Nra_ng=Nra-Nra_g;
else
    Nra_g=0;
    Nra_ng=Nra;
end

max_na=na;
Nra_a=Nra;
tmax = 5.4 * (10 ^-3);% TXOP duration
Vs=vs;an_st=vs;
f_subcarrier=78.125;% subcarrier spacing in kHz
n_tones =1000*(b/(Nru*f_subcarrier));% Number of tones
an_ap=ap_ant;

if(n_tones <=106) %  MU-MIMO is capable or not
    ap_ant=1;an_ap=ap_ant;
    csi_flag=0;
    flag_mimo=0;
end

time_csi=0;
t_ofdm=16;%microsec

[r]=find_r_ap(Nru,b,vs,mcs,t_ofdm,f_subcarrier); % Find Data-rate
[r_csi]=find_r_ap(74,b,vs,mcs,t_ofdm,f_subcarrier);% Fidn CSI data-rate
ysc=(b/f_subcarrier)*1000;
ysc_csi=(b/f_subcarrier)*1000;

r_ul=r;
pkt=PayLoad /r; %fraction of time
packet=PayLoad;
W_ul=W;
W_ng=W;W_g=W;
[cur_sim_na] = find_na_mu();% Find the actual Number fo packets aggregated

%% Other initialization
CWmin = W;
sim_time = 0;
pkts_drop_an_gm=[];pkts_drop_si_gm=[];drop_time_g_im=[];var_g_anam=[];total_delay_avgm=[];
pkts_drop_an_ngm=[];pkts_drop_si_ngm=[];drop_time_ng_im=[];var_ng_anam=[];total_delay_avngm=[];
per_sta_theory=[];per_sta_theory=[per_sta_theory 0];
throughput = [];bsrdelivery=[];bsr_simu=[];bsr_anal_mc=[];bsr_anal_fpa=[];
pkts_drop_an_ng=[];pkts_drop_si_ng=[];
throughputa_fpa = [];bsra=[];factor=[];index=1;throughputa_mc = [];pkts_drop_an_g=[];pkts_drop_si_g=[];
throughput=[throughput 0];throughputa_fpa=[throughputa_fpa 0];throughputa_mc=[throughputa_mc 0];
pkts_drop_an_g=[pkts_drop_an_g 0];
pkts_drop_si_g=[pkts_drop_si_g 0];
pkts_drop_an_ng=[pkts_drop_an_ng 0];
pkts_drop_si_ng=[pkts_drop_si_ng 0];
pkts_drop_an_gm=[pkts_drop_an_gm 0];delay_csi_g=0;delay_csi_ng=0;
pkts_drop_si_gm=[pkts_drop_si_gm 0];
drop_time_g_im=[drop_time_g_im 0];var_g_anam=[var_g_anam 0];
total_delay_avgm=[total_delay_avgm 0];

pkts_drop_an_ngm=[pkts_drop_an_ngm 0];
pkts_drop_si_ngm=[pkts_drop_si_ngm 0];
drop_time_ng_im=[drop_time_ng_im 0];var_ng_anam=[var_ng_anam 0];
total_delay_avngm=[total_delay_avngm 0];

%% Simulation
fprintf('Simulating...\n');   % kmin=10;

drop_time_g_is=[];drop_time_ng_is=[];tau_g_array=[];tau_ng_array=[];
delay_analyng111=[];drop_time_g_i=[];drop_time_ng_i=[];delay_analyngm=[];
theory_npkts=0;per_sta_packets_nodes=[];delay_analyg=[];delay_analyng=[];total_delay_n_ng=[];total_delay_avg=[];total_delay_avg_sim=[];total_delay_avg_sim_ng=[];total_delay_avg_ng=[];var_analy=[];var_g_ana=[];
tau_s=[];p_=[];tau_a=[];ps_ana=0;p_c_stas_the=[];p_c_stas_simu=[];total_delay_n_g=[];time_calc=0;total_stages_n=[];
tau_sg=[];tau_ag=[];delay_analyng=[delay_analyng 0];delay_analyngm=[delay_analyngm 0];
tau_sng=[];tau_ang=[];pc_sg=[];pc_sng=[];pc_ag=[];pc_ang=[];total_variance_ng=[];total_variance_g=[];
total_delay_avg_sim=[total_delay_avg_sim 0];total_delay_avg=[total_delay_avg 0];delay_analyng111=[delay_analyng111 0];
total_delay_avg_sim_ng=[total_delay_avg_sim_ng 0];total_delay_avg_ng=[total_delay_avg_ng 0];bsr_simu=[bsr_simu 0];var_analym=[];
bsr_anal_fpa=[bsr_anal_fpa 0];bsr_anal_mc=[bsr_anal_mc 0];total_variance_ng=[total_variance_ng 0];var_analy=[var_analy 0];var_analym=[var_analym 0];
var_g_ana=[var_g_ana 0];total_variance_g=[total_variance_g 0];
drop_time_g_i=[drop_time_g_i 0];
drop_time_ng_i=[drop_time_ng_i 0];
drop_time_g_is=[drop_time_g_is 0];
drop_time_ng_is=[drop_time_ng_is 0];
tau_g_array=[tau_g_array 0];
tau_ng_array=[tau_ng_array 0];energy_anal=[];energy_eff_anal=[];energy_simu=[];energy_eff_simu=[];

%% Control Frames Initialisation
time_calc22= times_exp_mu(Nru);
timeq=time_calc22;
timeq1=timeq;
time_csi=0;
t_tfr1= t_tfr;t_bsr1= t_bsr;t_trigger1= t_trigger;t_mu_data_ul1= t_mu_data_ul;t_msback1= t_msback;

%% Stations For loop
for num_stationb=k1
    drop_time_ng_isi=[];
    drop_time_g_isi=[];
    total_delay_array_g=[];total_delay_array_ng=[];total_sec_mom=[];bsr=[];

    tau_ni=[];p_i=[];tau_i=[];pc_i=[];through_simu = [];p_s_prob=[];n_pkts_stas=[];ddg=[];  dd_stages=0;  ddng=[];pc_ig=[];pc_ing=[];ddng_sec=[];
    tau_ig=[];tau_ing=[];time_csi=0;tot_del_all_ietr=[];tot_del_all_ietr_g=[];
    pkts_dis_inter12_g=[];pkts_dis_inter12_ng=[];

    twt_stas=[];
    if(twt_flag)% if TWT flag is enabled
        twt_stas=randi(T_max,1,num_stationb);% select the STAs based on t_Mean
    end

    for sta_iter = 1:1:n_iter
        n=num_stationb;n_act_stas=0;total_stas=n;
        if(twt_flag)
            n_th=(2*n)/(1+T_max);
        else
            n_th=n;
        end
        if(twt_flag)
            for ikl=1:1:n
                val=mod(sta_iter,twt_stas(ikl));
                if(val==0)
                    n_act_stas=n_act_stas+1;
                end
            end
            per_sta_packets=zeros(1,n_act_stas);nra=n;num_station=n;
        else
            n = num_stationb;per_sta_packets=zeros(1,num_stationb);nra=n;num_station=n;
        end

        if(flag_mimo)% Select the STAs capable of MU-MIMO
            vu=ceil(grouping_ratio*n); g=vu;
        else
            vu=0;g=0;
        end
        act_stas=n;
        
        if(n~=0)
            % CSI Calculation
            [time_csi] = get_csi_stage1(vu,r_csi,ysc_csi,an_ap,0);% Channel sounding Duration
            factor1=time_csi*lcsi*10e-6;
            factor(index)=1-factor1;
            if(factor(index)<0)
                factor(index)=0;
            end
            delay = zeros(1,num_station);stage= zeros(1,num_station);
            station_stage = zeros(1,num_station);
            next_tran_time = zeros(1,num_station);
            cw_stas = zeros(1,num_station);

            nCycles_noTx=0;
            nCycles_allSA=0;
            nCycles_t1=0;
            nCycles_allCollisions =0;
            nCycles_success_ra=0;
            group_flag_stas_all=zeros(nra,1);
            sim_time=0;

            suc_pkt = 0;
            collision_pkt = 0;
            cur_contending_stations = 1:n;
            n_contending_stations = zeros(tf_cycles, 1);
            n_transmitting_stations = zeros (tf_cycles,1);CW =zeros(n,1);OBO=zeros(n,1);total_delayg=[];total_delayng=[];
            total_stages=[];

            %% Initial back-off counter values initialisation
            for kkl=1:1:nra
                if(kkl<=g)
                    group_flag_stas_all(kkl,1)=1;
                    CW(kkl,1)=CWmin;
                    OBO(kkl,1) = floor(unifrnd (0,CWmin-1,1,1));
                else
                    CW(kkl,1)=CWmin;
                    OBO(kkl,1) = floor(unifrnd (0,CWmin-1,1,1));
                end
            end
            pkts_discarded_ng=0;n_suc_g=0;n_suc_ng=0;
            n_pkts_actual=[];pkts_discarded_g=0;n_suc=0;pkts_dis_inter1=0;
            drop_time_g=[];drop_time_ng=[];n_g_attempt=0;n_ng_attempt=0;
            for iteration = 1:tf_cycles
                ntx_g=0;ntx_ng=0;
                for kul=1:1:n
                    delay(1,kul)=delay(1,kul)+ timeq; % Delay parameter calculation
                end
                RU_status = zeros(Nra,lim);
                Group_status= zeros(Nra,lim);
                Coll_status= zeros(Nra,lim);
                cur_OBO = OBO(cur_contending_stations);

                %% Backoff Counter Decrementation
                for kkl=1:1:nra
                    if(group_flag_stas_all(kkl,1)==1)
                        cur_OBO(kkl,1) = cur_OBO(kkl,1) - Nra_g;
                    else
                        cur_OBO(kkl,1) = cur_OBO(kkl,1) - Nra_ng;
                    end
                end

                %% Finding the STAs that are successful
                OBO(cur_contending_stations) = cur_OBO;
                tx_nodes = find(OBO <= 0); %Set of nodes that transmit in this cycle
                n_tx = length(tx_nodes);
                for kkl=1:1:nra
                    if(group_flag_stas_all(kkl,1)==1)
                        if(cur_OBO(kkl,1) <=0)
                            n_g_attempt=n_g_attempt+1;
                        end
                    else
                        if(cur_OBO(kkl,1) <=0)
                            n_ng_attempt=n_ng_attempt+1;
                        end
                    end
                end
                nrag=g;nrang=nra-g;
                if(n_tx~=0)
                    for nodei=1:1:length(tx_nodes)
                        if(group_flag_stas_all(tx_nodes(nodei),1)==1)
                            ntx_g=ntx_g+1;
                        else
                            ntx_ng=ntx_ng+1;
                        end
                    end
                end
                n_transmitting_stations (iteration, 1) = n_tx;
                group_flag_stas=zeros(n_tx,1);


                ru_tx=[];nra_index=0;nra_array_group=[];
                for lll=1:1:Nra_g
                    nra_array_group=[nra_array_group lll];
                end
                nra_array_non_group=[];
                for lll=Nra_g+1:1:Nra
                    nra_array_non_group=[nra_array_non_group lll];
                end
                prob_g=g/nra;
                for kkl=1:1:n_tx
                    rand_n=rand;
                    if( group_flag_stas_all(tx_nodes(kkl),1)==1)
                        group_flag_stas(kkl,1)=1;
                    end
                end
                for kkk=1:1:n_tx
                    if group_flag_stas(kkk,1)==1 %Assigning rus to grouped nodes
                        value=randi([1,Nra_g],1,1);
                        nra_index=value;
                        ru_tx(kkk,1)=nra_array_group(nra_index);
                    else
                        if(Nra_ng~=0)
                            value=randi([1,Nra-Nra_g],1,1);
                            nra_index=value;
                            ru_tx(kkk,1)=nra_array_non_group(nra_index);
                        end

                    end
                end
                for i = 1:length(ru_tx)
                    RU_status(ru_tx(i,1),1) = RU_status(ru_tx(i,1),1)+ 1;
                    for j=2:1:lim
                        if (RU_status(ru_tx(i,1),j)==0)
                            RU_status(ru_tx(i,1),j) = tx_nodes(i,1);
                            if(group_flag_stas(i,1)==0)
                                Group_status(ru_tx(i,1),j)=1;%not grouped
                            else
                                Group_status(ru_tx(i,1),j)=2;%groupped
                            end
                            break;
                        end
                    end
                end
                ng_count=zeros(1,Nra);g_count=zeros(1,Nra);coll_flag=zeros(1,Nra);
                for i = 1:(Nra)
                    for j=2:1:lim
                        if( Group_status(i,j)~=0)
                            if(Group_status(i,j)==2)
                                g_count(1,i)=g_count(1,i)+1;
                            else
                                ng_count(1,i)=ng_count(1,i)+1;
                            end
                        else
                            break;
                        end
                    end
                end
                jj=1;
                for i = 1:(Nra)
                    if(g_count(1,i)>an_ap)%collision occured
                        if(Coll_status(i,jj)==0)
                            Coll_status(i,jj)=1;
                            jj=jj+1;
                        end
                    end
                    if(ng_count(1,i)>1)%collision occured
                        if(Coll_status(i,jj)==0)
                            Coll_status(i,jj)=2;
                            jj=jj+1;
                        end
                    end
                    if((g_count(1,i)~=0)&&(ng_count(1,i)~=0))%collision occured
                        if(Coll_status(i,jj)==0)
                            Coll_status(i,jj)=3;
                            jj=jj+1;
                        end
                    end
                end
                for i = 1:(Nra)
                    if(sum(Coll_status(i,:))>0)
                        coll_flag(1,i)=1;%collision occured on this RU
                    end
                end

                succ_stas=n_tx;col=0;
                for i = 1:length(ru_tx)
                    if(coll_flag(1,ru_tx(i))==1)
                        col=col+1;
                    end
                end

                succ_stas=n_tx-col;
                pc_i=[pc_i col/n_tx];
                p_s_prob=[p_s_prob succ_stas/nra];
                time_calc= times_exp_mu(((Nsa*an_ap)+succ_stas));

                %% Updating the Back-off stage and contention window parameters based on the STA success criteria
                successful_tx_nodes = [];nc_g=0;nc_ng=0;
                for i = 1:Nra
                    ru = i;
                    for tx_node_index=2:1:lim
                        node = RU_status(ru,tx_node_index);
                        if(node==0)
                            break;
                        end
                        if (coll_flag(1,ru) == 1) %a collision occurred on this RU, double CW
                            stage(1,node) = stage(1,node) +1;
                            if( stage(1,node)  <= m)
                                CW (node,1) = (2^stage(1,node))*W;
                            elseif( stage(1,node)  <= m_limit)
                                CW (node,1) = (2^m)*W;
                            else
                                if(group_flag_stas_all(node,1)==1)
                                    pkts_discarded_g=pkts_discarded_g+1;
                                    drop_time_g=[drop_time_g delay(1,node)];
                                else
                                    pkts_discarded_ng=pkts_discarded_ng+1;
                                    drop_time_ng=[drop_time_ng delay(1,node)];
                                end

                                stage(1,node) =0;
                                CW (node,1) = CWmin;
                                delay(1,node)=0;
                            end
                            if((group_flag_stas_all(node,1)==1))
                                nc_g=nc_g+1;
                            else
                                nc_ng=nc_ng+1;
                            end
                        else
                            successful_tx_nodes = [successful_tx_nodes node];
                            CW (node,1) = CWmin; %successful transmission occurred on this RU; reset CW
                            n_suc=n_suc +1;
                            stage(1,node) = 0;
                            if(group_flag_stas_all(node,1)==1)
                                [total_delayg]=[total_delayg delay(1,node)];
                                n_suc_g=n_suc_g +1;
                            else
                                [total_delayng]=[total_delayng delay(1,node)];
                                n_suc_ng=n_suc_ng +1;
                            end
                            [total_stages]=[total_stages stage(1,node)];
                            delay(1,node)=0;
                        end
                        OBO (node,1) = floor(unifrnd (0,CW (node,1)-1,1,1)); % Pick a new OBO
                    end
                end

                %******************************************************
                pkts_dummy=length(successful_tx_nodes);
                n_pkts_actual = [n_pkts_actual pkts_dummy];
                if(ntx_g~=0)
                    pc_ig=[pc_ig nc_g/ntx_g];
                end
                if(ntx_ng~=0)
                    pc_ing=[pc_ing nc_ng/ntx_ng];
                end
                [T1c] = times_exp_mu((pkts_dummy));

                through_simu =[through_simu (packet*((Nsa*an_ap)+pkts_dummy))];

                if(Nsa==0 && length(successful_tx_nodes)~=0)
                    nCycles_success_ra =nCycles_success_ra+1;
                end
            end
        end
if(~isempty(total_delayg))
        if(length(total_delayg)>0)
            ddg=[ddg mean(total_delayg)];
            tot_del_all_ietr_g=[tot_del_all_ietr_g total_delayg];
        end
else
    ddg=[ddg 0];tot_del_all_ietr_g=[tot_del_all_ietr_g 0];
end
        if(length(total_delayng)>0)
            ddng=[ddng mean(total_delayng)];
            tot_del_all_ietr=[tot_del_all_ietr total_delayng];
        end

        if(length(total_delayng)>0)
            asd=total_delayng .^3;
            ddng_sec=[ddng_sec mean(asd)];
        end

        dd_stages=[dd_stages mean(total_stages)];
        bsr_mean = mean(n_pkts_actual)/num_station;
        n_pkts_stas=[n_pkts_stas mean(n_pkts_actual)];

        bsr=[bsr bsr_mean];
    end



    %% Iteration End

    if(~isempty(drop_time_ng_isi))
        drop_time_ng_is=[drop_time_ng_is mean(drop_time_ng_isi)];
    else
        drop_time_ng_is=[drop_time_ng_is 0];

    end
    if(~isempty(drop_time_g_isi))
        drop_time_g_is=[drop_time_g_is mean(drop_time_g_isi)];
    else
        drop_time_g_is=[drop_time_g_is 0];

    end

    pkts_dis_simu_ng=mean(pkts_dis_inter12_ng);
    pkts_dis_simu_g=mean(pkts_dis_inter12_g);

    per_sta_packets=per_sta_packets/n_iter;
    per_sta_packets=per_sta_packets/num_station;
    total_delay_array_g=[total_delay_array_g ddg];
    total_delay_array_ng=[total_delay_array_ng ddng];
    total_sec_mom=[total_sec_mom ddng_sec];
    mean_ng=mean(total_delay_array_ng);
    var_ng=0;co=1;
    for hj=1:1:length(total_delay_array_ng)
        var_ng= var_ng+(total_delay_array_ng(hj)-mean_ng)^2;
        co=co+1;
    end
    var_ng=sqrt(var_ng/(co-1));
    var_direct_formula=sqrt(var(tot_del_all_ietr));
    if(~isempty(tot_del_all_ietr))
        var_ng= sqrt(var(tot_del_all_ietr));
    else
        var_ng=0;
    end
    Simu_varr=var_ng;
    if(~isempty(tot_del_all_ietr_g))
        var_g= sqrt(var(tot_del_all_ietr_g));
        total_variance_g=[total_variance_g var_g];
    else
        var_g=0;  total_variance_g=[total_variance_g var_g];
    end
    if(~isempty( ddg))
        total_delay_n_g=[total_delay_n_g mean(ddg)];
    end
    if(~isempty( ddng))
        total_delay_n_ng=[total_delay_n_ng mean(ddng)];
    end
    total_variance_ng=[total_variance_ng var_ng];
    total_stages_n=[total_stages_n mean(dd_stages)];
    if(n~=0 )

        %% Throughput and RA Capability Calculation
        if(csi_flag)
            bsr_simu=[bsr_simu factor(index)*mean(bsr)];
            throughput = [throughput factor(index)*(packet*cur_na*((Nsa*an_ap)+((mean(n_pkts_stas)))))/(timeq)];
            [th,ra]= func_new_fpa(m_limit,n_th,Nsa*an_ap); a_ng=gamma_nga/Nra_ng;
            if(grouping_ratio==0)
                ra=(nchoosek(nra,1))*(a_ng)*((1-a_ng)^(nra-1))*Nra_ng;
            else
                ra=ra*factor(index);
            end
            tau_g_array=[tau_g_array tau_ga];
            tau_ng_array=[tau_ng_array tau_nga];
            throughputa_fpa = [throughputa_fpa, factor(index)*cur_na*th];
            bsr_anal_fpa=[bsr_anal_fpa ra];
            [th1,ra1]= func_new_mc(n);
            a_ng=gamma_nga/Nra_ng;
            if(grouping_ratio==0)
                ra=(nchoosek(nra,1))*(a_ng)*((1-a_ng)^(nra-1))*Nra_ng;
            else
                ra=ra*factor(index);
            end
            bsr_anal_mc=[bsr_anal_mc ra1];
            throughputa_mc = [throughputa_mc, factor(index)*cur_na*th1];                %
        else
            bsr_simu=[bsr_simu mean(bsr)];
            n_suc=mean(n_pkts_stas);
            throughput = [throughput (packet*cur_na*((Nsa*an_ap)+((mean(n_pkts_stas)))))/(timeq)];
            [th,ra]= func_new_fpa(m_limit,n,Nsa*an_ap);
            tau_g_array=[tau_g_array tau_ga];
            tau_ng_array=[tau_ng_array tau_nga];
            throughputa_fpa = [throughputa_fpa, cur_na*th];
            bsr_anal_fpa=[bsr_anal_fpa ra];
            [th1,ra1]= func_new_mc(n);bsr_anal_mc=[bsr_anal_mc ra1];
            throughputa_mc = [throughputa_mc, cur_na*th1];
        end

        %% Delay Moments Calculation
        if(Nra_g~=0)
            if(gamma_ga~=0)
                delayxng1112=find_moments_mean_alt(W,m_limit,m,gamma_ga,Nra_g);
                delayxg=delayxng1112;
            else
                std_another=0;delayxg=0;
            end
            total_delay_avg=[total_delay_avg delayxg];delayxgm=0;
            if(gamma_g_mc~=0)
                delayxng1112m=find_moments_mean_alt(W,m_limit,m,gamma_g_mc,Nra_g);
                delayxgm=delayxng1112m;
            end
            %%% find delay with csi
            if(csi_flag)
            delay_csi_g=find_tot_delay(delayxgm*(10^-6),lcsi,time_csi*(10^-6));
            
            else
               delay_csi_g=0;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
total_delay_avgm=[total_delay_avgm delay_csi_g];
        delayxng111=find_moments_mean_alt(W,m_limit,m,gamma_nga,Nra_ng);

        if(num_station~=g)
            delayxng=delayxng111;
        else
            delayxng=0;
        end
        if(csi_flag)
            delay_csi_ng=find_tot_delay(delayxng*(10^-6),lcsi,time_csi*(10^-6));

        else
            delay_csi_ng=0;
        end
        delay_analyng=[delay_analyng delay_csi_ng];
        gamma_ngam=gamma_ng_mc;
        delayxng111m=find_moments_mean_alt(W,m_limit,m,gamma_ngam,Nra_ng);%mean(pc_ing)
        if(num_station~=g)
            delayxngm=delayxng111m;
        else
            delayxngm=0
        end
        delay_analyngm=[delay_analyngm delayxngm];

       %% Energy Calculation
        n_att_groupes=(n_g_attempt+n_ng_attempt)/tf_cycles;
        if(csi_flag)
            [eff,energy_e]=find_energy_csi(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,t_g,Nsa*an_ap,lcsi,time_csi,n_suc);
        else
            [eff,energy_e]=find_energy(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,t_g,Nsa*an_ap,n_suc);
        end
        energy_anal=[energy_anal energy_e];
        energy_eff_anal=[energy_eff_anal eff];



        %% Energy effeciency
        if(csi_flag)
            [effsim,energy_esim]=find_energy_csi(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,n_att_groupes,Nsa*an_ap,lcsi,time_csi,n_suc);
        else
            [effsim,energy_esim]=find_energy(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,n_att_groupes,Nsa*an_ap,n_suc);
        end
        energy_simu=[energy_simu energy_esim];
        energy_eff_simu=[energy_eff_simu effsim];
    end
    index=index+1;

    %% For printing purpose

    num_station

    if(~isempty( total_delay_n_g))
        Delay_Simg= total_delay_n_g(length(total_delay_n_g));
        s1_g= total_delay_n_g(length(total_delay_n_g));
    else
        Delay_Simg = 0;
        s1_g=0;
    end
    if(~isempty( delay_analyg))
        Delay_Anag= delay_analyg(length(delay_analyg))
    end

    if(~isempty(total_delay_n_ng))
        Delay_Simng= total_delay_n_ng(length(total_delay_n_ng));
    end
    if(~isempty(delay_analyng))
        Delay_Anang= delay_analyng(length(delay_analyng));
    end
    if(~isempty(delay_analyng111))
        delay_other=delay_analyng111(length(delay_analyng111));
    end
    s1_ng=0;
    if(~isempty(total_delay_n_ng))
        s1_ng=total_delay_n_ng(length(total_delay_n_ng));
    end
    if(Nra_ng~=0)
        s1_avg=(s1_g);
    else
        s1_avg=(s1_g);
    end
    del_g_simu=find_tot_delay(s1_avg*(10^-6),lcsi,time_csi*(10^-6));
    if(~isempty(s1_avg))
        if(csi_flag)
        total_delay_avg_sim=[total_delay_avg_sim (del_g_simu)];%+((s1_avg)*(10^-6))];
        else
              total_delay_avg_sim=[total_delay_avg_sim ((s1_avg)*(10^-6))];
        end
    end
 
    del_ng_simu=find_tot_delay(s1_ng*(10^-6),lcsi,time_csi*(10^-6));
    if(~isempty(s1_ng))
        if(csi_flag)
            total_delay_avg_sim_ng=[total_delay_avg_sim_ng (del_ng_simu)];%+((s1_ng)*(10^-6))];
        else
            total_delay_avg_sim_ng=[total_delay_avg_sim_ng ((s1_ng)*(10^-6))];
        end
    end
end

total_delay_avngm = total_delay_avngm*(10^-6);
total_delay_avg=total_delay_avg*(10^-6);
total_delay_avg_ng=total_delay_avg_ng*(10^-6);
delay_analyng111=delay_analyng111*(10^-6);
delay_analyngm=delay_analyngm*(10^-6);
throughputa_fpa_2=[];throughputa_mc_2=[];
throughputa_fpa_2=[throughputa_fpa_2 throughputa_fpa(1)];
throughputa_mc_2=[throughputa_mc_2 throughputa_mc(1)];

%% Output Metrics
Throughput_Simu=throughput
Throughput_Ana=throughputa_mc
RA_Capability_Simu=bsr_simu
RA_Capability_Ana=bsr_anal_mc
if(csi_flag)
Delay_Simu=total_delay_avg_sim
Delay_Ana=total_delay_avgm
end
if(csi_flag==0)
Delay_Ana=delay_analyngm
Delay_Simu=total_delay_avg_sim_ng
end
PowerSimu=energy_simu
Power_Ana=energy_anal
Effe_Simu=energy_eff_simu
Effe_Ana=energy_eff_anal


