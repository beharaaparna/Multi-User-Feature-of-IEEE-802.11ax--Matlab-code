function [th,na_out]=find_dyn_thr_poisson_edited_cc(Nra,packet,timeq,lambda,m_limit,an_ap,n,m)
global lambda;global na_max;global new_stas;global dyn_Nra;global dyn_flag;global buffer_max;global flag_finite_buffer;global max_sched_slots;
if(dyn_flag)
    Mx=700;
else
    Mx=1;
end
Nra_dyn=Nra*ones(1,Mx+11111113);
Nsa_dyn=zeros(1,Mx+111113);
New_STAs=n*ones(1,Mx+111111113);
thr_array=[];slots_sch_array=[];
n_pkts_next_sch=zeros(1,Mx+111111131);
p_suc=zeros(1,Mx+1111131);
    na_theory=[];
    
for tf=1:1:Mx
    [th,n_suc,n_pkts_suc,n_pkts_trans,gamma]= func_new_fpa_nonsat_dyn_cc(m_limit,New_STAs(1,tf),Nsa_dyn(1,tf),Nra_dyn(1,tf));
    thr_array=[thr_array th];
 
    if(~dyn_flag)
        Nsa_dyn(1,tf+1)=0;Nra_dyn(1,tf+1)=Nra;New_STAs(1,tf+1)=n;
         pl=max(n_pkts_suc,1);
%            if(flag_finite_buffer)
%           na_fresh=min(buffer_max,pl);
%            else
%                na_fresh=pl;
%            end
             na_theory=[na_theory n_pkts_trans];%min(na_max,na_fresh)];

    else
        pl=max(n_pkts_suc,1);
        if(flag_finite_buffer)
            na_fresh=min(buffer_max,pl);

            na_left=min(buffer_max,na_fresh+na_fresh-min(na_max,na_fresh));
        else
            na_left=n_pkts_suc;%-n_pkts_trans;
        end
        left_pkts_tracker=na_left;slots_Sched1=0;
        if(left_pkts_tracker<=0)
            na_theory=[na_theory n_pkts_trans];%min(na_max,na_left)];
        end
        while(left_pkts_tracker>0)
            slots_Sched1=slots_Sched1+1;
            pll=min(left_pkts_tracker,na_max);%min(na_max,left_pkts_tracker);
            left_pkts_tracker=left_pkts_tracker-pll;
            na_theory=[na_theory pll];
%              if(slots_Sched1>=max_sched_slots)
%                 break;
%              end
        end
         total_slots_sch=slots_Sched1;
%         total_slots_sch=min(max_sched_slots,slots_Sched1);
        slots_Sched=slots_Sched1;
        
%         p_suc(1,tf+total_slots_sch)=n_suc/New_STAs(1,tf);
        slots_sch_array=[slots_sch_array slots_Sched];

         x=tf+1;
        while(slots_Sched>1)
            if(n_suc>100)
                here=1;
            end
            x
            Nsa_dyn(1,x)
            n_suc
            Nsa_dyn(1,x)=Nsa_dyn(1,x)+n_suc;
            if(Nsa_dyn(1,x)>Nra)
                Nsa_dyn(1,x)=Nra;
            end
            New_STAs(1,x)=New_STAs(1,x)-n_suc;
            Nra_dyn(1,x)=Nra-Nsa_dyn(1,x);
            x=x+1;
            slots_Sched=slots_Sched-1;
%             n_pkts_next_sch(1,tf+total_slots_sch)=n_pkts_next_sch(1,tf+total_slots_sch)+find_n_pkts_accumulated2(lambda,timeq*(10^-6));
        end

%         Nsa_dyn(1,tf+1)=0;
%         Nra_dyn(1,tf+1)=Nra;
%         New_STAs(1,tf+1)=n;
    end
    if(tf==15)
        break;
    end
end
th=mean(thr_array);
na_out=mean(na_theory);