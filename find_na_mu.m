function [na] = find_na_mu()
global max_na;
global tmax;
global cur_na;
global r_ul;
global mh;global md;global PayLoad;
md=32;%in bits
mh=320;%in bits;
na_temp = floor((r_ul * tmax * (10 ^6))/(mh+md+PayLoad));
if (na_temp <= max_na)
    na = na_temp;
else
    na = max_na;
    while(1)
        if((na/(r_ul*(10 ^6)))<= tmax)
            break;
        else
            na=na-1;
            if(na==1)
                break;
            end
        end
    end
end
cur_na=na;
end