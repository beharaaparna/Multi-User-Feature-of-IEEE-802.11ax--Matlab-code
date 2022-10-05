clear all;
clc;
close all;

%% Global Parameters Declaration 
global r;
global act_stas;
global gamma_ga;global gamma_nga;
global PayLoad;global na_max;global K_scale;global alpha;
global W_ul;global W_g;global W_ng;
global Nra_a;global Nsa;global b;
global nra;global Nra_ng;
global m;global an_ap;global Nru;global g;global Nra_g;global an_ap;global Vs;global an_st;global max_na;
global mcs;global cur_na;
global tmax;global m_ul;global r_ul;global Nra;global m_limit;global timeq;
global gamma_g_mc;
global gamma_ng_mc;
global tau_nga;global tau_ga;
global t_g;global fn_array;
global t_tfr;global t_bsr;global t_trigger;global t_mu_data_ul;global t_msback;
global t_tfr1;global t_bsr1;global t_trigger1;global t_mu_data_ul1;global t_msback1;
global lambda_max;
global afille;global K_scale;global alpha;global lambda;

%% Simulation parameters as given in Table.1
csi_flag=0;% Flag to enable Channel Sounding
lambda_max=80;% Packet Arrival Rate
prob_same_tid=1;%Pobability that the packets are of the same traffic identifier
twt_flag=0;% Flag to enable TWT
lambda=lambda_max;
[K_scale,alpha]= find_pareto_lambda(lambda);% Find Paretto Parameters
T_max=15;% TWT t_mean
flag_mimo=0; % Flag to enable MU-MIMO
lcsi=5;% Channel Sounding Interval
ap_ant=1; % Number fo AP antennas
Nru=16;% Number of RUs
mcs=11;% Modulation & Coding Index
Nra=Nru;% Number of RA RUs
Nsa=Nru-Nra;% Number of SA RUs
grouping_ratio=0.01;% Percentage of STAs dedicated for MU-MIMO
afille=fopen("thr_test.txt",'w');% File to write the inbtermediate results
tf_cycles =10000;% Number of time-slots (TF-Cycles) per iteration considered
n_iter =10;% Number of Iterations
k1=[50 100:100:2000];
m_limit=5;% Maximum Retry Limit
PayLoad = 12000;%paylaod in bits
W=64;% Minimum Contention Window Size
na=256;% Maximum number fo packets that can be aggregated
m=5;
m_ul=m;
vs=1;% Number of antennas at the STA
b=160;% Bandwidth in MHz

if(flag_mimo) % MIMO Flag if enabled, selects the RUs for MU-MIMO transmissions
    Nra_g=ceil(Nra*grouping_ratio);
    Nra_ng=Nra-Nra_g;
else
    Nra_g=0;
    Nra_ng=Nra;
end

max_na=na;
Nra_a=Nra;
tmax = 5.4 * (10 ^-3);% TXOP duration
Vs=vs;an_st=vs;
lim=10000;
f_subcarrier=78.125;% subcarreir spacing in Khz
n_tones =1000*(b/(Nru*f_subcarrier));% Calculating the number of tones
an_ap=ap_ant;
if(n_tones <=106)
    ap_ant=1;an_ap=ap_ant;% If number of tones are less than 106, MU-MIMO transmissiosn not feasible
end

time_csi=0;
t_ofdm=16;%OFDM Symbol duration microsec
[r]=find_r_ap(Nru,b,vs,mcs,t_ofdm,f_subcarrier);% Data-rate
[r_csi]=find_r_ap(74,b,vs,mcs,t_ofdm,f_subcarrier);% Data-Rate fro CSI transmissions
ysc=(b/f_subcarrier)*1000;
ysc_csi=(b/f_subcarrier)*1000;
r_ul=r;
pkt=PayLoad /r; %fraction of time
packet=PayLoad;
W_ul=W;
W_ng=W;W_g=W;
[cur_sim_na] = find_na_mu();% Calculate the actual number fo packets aggregated
na_max=cur_sim_na;

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
pkts_drop_an_gm=[pkts_drop_an_gm 0];
pkts_drop_si_gm=[pkts_drop_si_gm 0];
drop_time_g_im=[drop_time_g_im 0];var_g_anam=[var_g_anam 0];
total_delay_avgm=[total_delay_avgm 0];
pkts_drop_an_ngm=[pkts_drop_an_ngm 0];
pkts_drop_si_ngm=[pkts_drop_si_ngm 0];
drop_time_ng_im=[drop_time_ng_im 0];var_ng_anam=[var_ng_anam 0];
total_delay_avngm=[total_delay_avngm 0];

%% Simulation
fprintf('Simulating...\n');   

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
tau_ng_array=[tau_ng_array 0];energy_anal=[0];energy_eff_anal=[0];energy_simu=[0];energy_eff_simu=[0];

%% Control Frame Initialisation
time_calc22= times_exp_mu(Nru);
timeq=time_calc22;
timeq1=timeq;
t_tfr1= t_tfr;t_bsr1= t_bsr;t_trigger1= t_trigger;t_mu_data_ul1= t_mu_data_ul;t_msback1= t_msback;

%% For loop for number fo STAs
for num_stationb=k1
    drop_time_ng_isi=[];
    drop_time_g_isi=[];
    total_delay_array_g=[];total_delay_array_ng=[];total_sec_mom=[];bsr=[];
    tau_ni=[];p_i=[];tau_i=[];pc_i=[];through_simu = [];p_s_prob=[];n_pkts_stas=[];ddg=[];  dd_stages=0;  ddng=[];pc_ig=[];pc_ing=[];ddng_sec=[];
    tau_ig=[];tau_ing=[];time_csi=0;tot_del_all_ietr=[];tot_del_all_ietr_g=[];
    pkts_dis_inter12_g=[];pkts_dis_inter12_ng=[];

    twt_stas=[];
    if(twt_flag)
        twt_stas=randi(T_max,1,num_stationb);% Based on TWT flag,selects the active STAs for the current time-slot
    end

    %% Iterations starts 
    for sta_iter = 1:1:n_iter
        n=num_stationb;n_act_stas=0;total_stas=n;
        %% TWT
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
            n = n_act_stas
            per_sta_packets=zeros(1,n_act_stas);num_station=n;
        else
            n = num_stationb;per_sta_packets=zeros(1,num_stationb);num_station=n;
        end

        %% MU-MIMO
        if(flag_mimo)
            vu=ceil(grouping_ratio*n); g=vu;
        else
            vu=0;g=0;
        end
        act_stas=n;% Actual STAs woke up in the current time-slot

        if(n~=0)
            %% CSI Calculation
            [time_csi] = get_csi_stage1(vu,r_csi,ysc_csi,an_ap,0);% CSI duration
            factor1=time_csi*lcsi*10e-6;
            factor(index)=1-factor1;% CSI overhead
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
            group_flag_stas_all=zeros(n,1);
            sim_time=0;
            suc_pkt = 0;
            collision_pkt = 0;
            cur_contending_stations = 1:n;
            n_contending_stations = zeros(tf_cycles, 1);
            n_transmitting_stations = zeros (tf_cycles,1);CW =zeros(n,1);OBO=zeros(n,1);total_delayg=[];total_delayng=[];
            total_stages=[];total_delayg_drop=[];

            %% Initial Back-Off Counter Initialisation
            for kkl=1:1:n
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
            drop_time_g=[];drop_time_ng=[];n_g_attempt=0;n_ng_attempt=0; n_pkts_buff=zeros(1,n);
            avg_pkt_buffer=[]; pkts_trans=[];len_cont_stas=[];lambda_array=[];

            for kl=1:1:n
                lambda_array=[lambda_array lambda];
            end

            %% Time-slots (TF-Cycles) advancements
            for iteration = 1:tf_cycles
                lam=[];
                if(iteration  ~= 0)
                    cur_contending_stations=[];lam=[];
                    for ik=1:1:n
                        [K,ALP]=find_pareto_lambda(lambda_array(ik));% Generate the packets using Pareto distributions
                        lam = [lam generate_pareto_rnd(K,ALP,1)];
                    end

                    for qq=1:1:n
                        if(lam(qq) < (timeq*(10^(-6))))
                            n_pkts_buff(qq)= n_pkts_buff(qq)+1; % Based on Paretto-based packet arrvial rate, find the number of packets present in each STA buffer
                        end
                    end

                    nra=0;
                    for qq=1:1:n
                        if(n_pkts_buff(qq)>0)
                            nra=nra+1;% Based on teh buffer status, find teh number of contending STAs for the current time-slot
                            cur_contending_stations = [cur_contending_stations qq];
                        end
                    end
                else
                    nra=n;
                end
                len_cont_stas=[len_cont_stas length(cur_contending_stations)];
                avg_pkt_buffer=[nanmean(n_pkts_buff)];
                ntx_g=0;ntx_ng=0;

                for kul=1:1:nra
                    noo=cur_contending_stations(kul);
                    delay(1,noo)=delay(1,noo)+ timeq;%% Delay variable
                end

                %% RU arrays used to store the RU status
                RU_status = zeros(Nra,lim);
                Group_status= zeros(Nra,lim);
                Coll_status= zeros(Nra,lim);
                cur_OBO = OBO(cur_contending_stations);

               %% Back-off Decrementation process
                for kkl=1:1:nra
                    if(group_flag_stas_all(kkl,1)==1)
                        cur_OBO(kkl,1) = cur_OBO(kkl,1) - Nra_g;
                    else
                        cur_OBO(kkl,1) = cur_OBO(kkl,1) - Nra_ng;
                    end
                end
                %% Finding the Succesful STAs for teh current TF-Cycle
                OBO(cur_contending_stations) = cur_OBO;
                tx_nodes = find(OBO(cur_contending_stations) <= 0); %Set of nodes that transmit in this cycle
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
%% Updating the back-off stage and contention window based on the success/failure criteria of each STA
                successful_tx_nodes = [];nc_g=0;nc_ng=0;disc_nodes=[];
                for i = 1:Nra
                    ru = i;
                    for tx_node_index=2:1:lim
                        node = RU_status(ru,tx_node_index);
                        if(node==0)
                            break;
                        end
                        nn_node=cur_contending_stations(node);
                        if (coll_flag(1,ru) == 1) %a collision occurred on this RU, double CW
                            stage(1,nn_node) = stage(1,nn_node) +1;
                            if( stage(1,nn_node)  <= m)
                                CW (nn_node,1) = (2^stage(1,nn_node))*W;
                            elseif( stage(1,nn_node)  <= m_limit)
                                CW (nn_node,1) = (2^m)*W;
                            else
                                if(group_flag_stas_all(nn_node,1)==1)
                                    pkts_discarded_g=pkts_discarded_g+1;
                                    drop_time_g=[drop_time_g delay(1,nn_node)];
                                else
                                    pkts_discarded_ng=pkts_discarded_ng+1;
                                    drop_time_ng=[drop_time_ng delay(1,nn_node)];
                                end
                                [total_delayg_drop]=[total_delayg_drop delay(1,nn_node)];
                                disc_nodes=[disc_nodes nn_node];
                                stage(1,nn_node) =0;
                                CW (nn_node,1) = CWmin;
                                delay(1,nn_node)=0;
                            end
                            if((group_flag_stas_all(nn_node,1)==1))
                                nc_g=nc_g+1;
                            else
                                nc_ng=nc_ng+1;
                            end
                        else
                            successful_tx_nodes = [successful_tx_nodes nn_node];
                            CW (nn_node,1) = CWmin; %successful transmission occurred on this RU; reset CW
                            n_suc=n_suc +1;
                            stage(1,nn_node) = 0;
                            if(group_flag_stas_all(nn_node,1)==1)
                                [total_delayg]=[total_delayg delay(1,nn_node)];
                                n_suc_g=n_suc_g +1;
                            else
                                [total_delayng]=[total_delayng delay(1,nn_node)];
                                n_suc_ng=n_suc_ng +1;
                            end
                            [total_stages]=[total_stages stage(1,nn_node)];
                            delay(1,nn_node)=0;
                        end
                        OBO (nn_node,1) = floor(unifrnd (0,CW (nn_node,1)-1,1,1)); % Pick a new OBO
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

                %% Successful transmissions
                for kkl=1:1:length(successful_tx_nodes)
                    nodee=successful_tx_nodes(kkl);
                    name_name=nodee;

                    if(rand<prob_same_tid)
                        vv=n_pkts_buff(name_name);
                        pkts_trans=[pkts_trans min(na_max,vv)];
                        if(vv <= na_max)
                            n_pkts_buff(name_name)= 0;
                        else
                            n_pkts_buff(name_name)= 0;% Update the buffer packet status of each STA
                        end
                    else
                        vv=n_pkts_buff(name_name);
                        pkts_trans=[pkts_trans 1];
                        n_pkts_buff(name_name)= 0;
                    end
                end

                %% Discarded Packets
                for kkl=1:1:length(disc_nodes)%disc_nodes
                    nodee=disc_nodes(kkl);
                    name_name=nodee;
                    n_pkts_buff(name_name)=0;
                end
                if(isempty(pkts_trans))
                    na_curr=0;
                else
                    gg=ceil(length(pkts_trans)/2);
                    na_curr=nanmean(pkts_trans(gg:length(pkts_trans)));
                end
                thr=(packet*((Nsa*an_ap)+(na_curr*(mean(n_pkts_actual)))))/(timeq);
                %% Logging the intermediate Outputs
                afille=fopen("thr_test.txt",'a');
                fprintf(afille,"TF_Cycles: %d, Num_Contend_STAs: %d, Throughput: %f \n",iteration, length(cur_contending_stations), thr );
                fclose(afille);
                fprintf("TF_Cycles: %d, Num_Contend_STAs: %d, Throughput: %f, Avg_Packet_length: %f \n",iteration, length(cur_contending_stations), thr, na_curr);

            end % end tf cycle
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
%% End of teh Iteration

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
        %% Throughput Calculation
        if(csi_flag)
            bsr_simu=[bsr_simu factor(index)*mean(bsr)];
            throughput = [throughput factor(index)*(packet*cur_na*((Nsa*an_ap)+((mean(n_pkts_stas)))))/(timeq)];
            [th,ra]= func_new_fpa(m_limit,n,Nsa*an_ap); a_ng=gamma_nga/Nra_ng;
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
            thr_simu= (packet*((Nsa*an_ap)+((mean(n_pkts_stas)))))/(timeq);
            gg=ceil(length(pkts_trans)/2);
            throughput = [throughput nanmean(pkts_trans(gg:length(pkts_trans)))*thr_simu];
            [th,ra,delay1,gamma,pkts_txted]= func_new_fpa_nonsat_pareto(m_limit,n,Nsa*an_ap);
            tau_g_array=[tau_g_array tau_ga];
            tau_ng_array=[tau_ng_array tau_nga];
            na_theory=pkts_txted;
            gg=ceil(length(pkts_trans)/2);
            na_simu= nanmean(pkts_trans(gg:length(pkts_trans)));
            throughputa_fpa = [throughputa_fpa, na_theory*th];
            bsr_anal_fpa=[bsr_anal_fpa ra];
        end
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
            total_delay_avgm=[total_delay_avgm delayxgm];
        end

        delayxng111=find_moments_mean_alt(W,m_limit,m,gamma_nga,Nra_ng); % Delay Moments

        if(num_station~=g)
            delayxng=delayxng111;
        else
            delayxng=0;
        end

        delay_analyng=[delay_analyng delayxng];
        gamma_ngam=gamma_ng_mc;
        delayxng111m=find_moments_mean_alt(W,m_limit,m,gamma_ngam,Nra_ng);%mean(pc_ing)
        if(num_station~=g)
            delayxngm=delayxng111m;
        else
            delayxngm=0;
        end
        delay_analyngm=[delay_analyngm delayxngm];

      %% Energy Calculation
        n_att_groupes=(n_g_attempt+n_ng_attempt)/tf_cycles;
        if(csi_flag)
            [eff,energy_e]=find_energy_csi(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,t_g,Nsa*an_ap,lcsi,time_csi);
        else
            eff=0;energy_e=0;
        end
        energy_anal=[energy_anal energy_e];
        energy_eff_anal=[energy_eff_anal eff];
      %% Energy effeciency
        if(csi_flag)
            [effsim,energy_esim]=find_energy_csi(total_stas,twt_flag,nra,throughputa_fpa(length(throughputa_fpa)),timeq,1,n_att_groupes,Nsa*an_ap,lcsi,time_csi);
        else
            effsim=0;energy_esim=0;
        end
        energy_simu=[energy_simu energy_esim];
        energy_eff_simu=[energy_eff_simu effsim];
    end
    index=index+1;

    %For printing purpose

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
    if(~isempty(s1_avg))
        total_delay_avg_sim=[total_delay_avg_sim s1_avg];
    end

    if(~isempty(s1_ng))
        total_delay_avg_sim_ng=[total_delay_avg_sim_ng s1_ng];
    end

end

Thr_Simulation=throughput
Thr_Analytical=throughputa_fpa

total_delay_avgm = total_delay_avgm*(10^-6);
total_delay_avngm = total_delay_avngm*(10^-6);
total_delay_avg=total_delay_avg*(10^-6);
total_delay_avg_sim=total_delay_avg_sim*(10^-6);
total_delay_avg_ng=total_delay_avg_ng*(10^-6);
total_delay_avg_sim_ng=total_delay_avg_sim_ng*(10^-6);
delay_analyng=delay_analyng*(10^-6);
delay_analyng111=delay_analyng111*(10^-6);
delay_analyngm=delay_analyngm*(10^-6);
throughputa_fpa_2=[];throughputa_mc_2=[];
throughputa_fpa_2=[throughputa_fpa_2 throughputa_fpa(1)];
throughputa_mc_2=[throughputa_mc_2 throughputa_mc(1)];
