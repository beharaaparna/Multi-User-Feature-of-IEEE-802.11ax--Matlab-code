function [r]=find_r_ap(nra,b,vs,mcs,t_symbol,f_spacing)
r=0;
n_tones=b/nra;
n_tones=n_tones/f_spacing;
n_tones=n_tones*1000;
if(mcs==6)
     r=6*n_tones*(3/4)*(vs/t_symbol);
elseif (mcs ==11)
    r=10*n_tones*(5/6)*(vs/t_symbol);
     elseif (mcs ==0)
    r=2*n_tones*(1/4)*(vs/t_symbol);
    elseif (mcs ==1)
    r=2*n_tones*(1/2)*(vs/t_symbol);
    elseif (mcs ==2)
    r=2*n_tones*(3/4)*(vs/t_symbol);
    elseif (mcs ==3)
    r=4*n_tones*(1/4)*(vs/t_symbol);
    elseif (mcs ==4)
    r=4*n_tones*(3/8)*(vs/t_symbol);
    elseif (mcs ==5)
    r=6*n_tones*(2/3)*(vs/t_symbol);
    elseif (mcs ==7)
    r=6*n_tones*(5/6)*(vs/t_symbol);
    elseif (mcs ==8)
    r=8*n_tones*(3/4)*(vs/t_symbol);
    elseif (mcs ==9)
    r=8*n_tones*(5/6)*(vs/t_symbol);
else
    r=10*n_tones*(3/4)*(vs/t_symbol);
  
end
    