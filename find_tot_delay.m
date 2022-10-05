function delay_csi=find_tot_delay(delay,lcsi,time_csi)
flag=0;
for k=1:1:1000000
    first=(k-1)*(1/lcsi);
    final=k*(1/lcsi);
    vv=(k*time_csi)+delay;
    if((vv > first) && (vv < final))
        flag=1;break;
    else
        here=1;
    end
end
if(flag==1)
delay_csi=((k)*time_csi) + delay;
end
delay_csi=delay_csi*(10^6);
end