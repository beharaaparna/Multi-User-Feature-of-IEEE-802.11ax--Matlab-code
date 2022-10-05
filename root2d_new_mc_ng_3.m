function F = root2d_new_mc_ng_3(x)
global W_ng;
global W_g;
global Nra_ng;
global m_ul;
global an_ap;
global g;
global nra;
global m_limit;

sum=0;sum12=0;
Nra=Nra_ng;w=W_ng;r=Nra;mbar=m_ul;m=m_limit;q=x(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b_1=(3*r)-2;
% b_2=(1-(q^(m+1)));
% b_3=1-q;
% if(b_2~=0)
%     b=b_1*b_2/b_3;
% else
%     b=0;
% end
% 
% sum1=0;
% rr=2*x(1);
% for i=0:1:mbar
%     if rr==0
%         sum1=sum1+0;
%     else
%         sum1=sum1+((rr)^i);
%     end
% end
% a1=w*sum1;
% 
% 
% sum2=0;
% for i=mbar+1:1:m
%     sum2=sum2+((x(1))^i);
% end
% a2=w*sum2*(2^mbar);
% 
% a=a1+a2;
% c=a;
f=m_limit-mbar;
sum1=0;
for i=0:1:m
    sum1=sum1+(q^i);
end

sum12=0;
for i=0:1:m
    if(i<mbar)
        w_eff=(2^i)*w;
    else
        w_eff=(2^mbar)*w;
    end
    term=(q^i)/w_eff;
    sum12=sum12+term;
end


sum123=0;
for i=0:1:m
    if(i<mbar)
        w_eff=(2^i)*w;
    else
        w_eff=(2^mbar)*w;
    end
    term=w_eff*(q^i);
    sum123=sum123+term;
end
sum123=sum123/(2*r);

sum4=0;
for i=0:1:m
    sum4=sum4+(q^i);
end
sum4=(r-2)*sum4;
sum4=sum4/(2*r);
sum=sum1+sum12+sum123+sum4;

b_oo=1/sum;

if(x(1)~=1)
    tau=b_oo *(1-(x(1)^(m+1)))/(1-x(1));
else
    tau=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_ng=tau/Nra;
sum=0;sum12=0;
Nra=Nra_ng;

if(nra-g>1)
    c=((1-a_ng)^(nra-g-1));
else
    c=1;
end

p1=c;

p=1-p1;
F(1)= x(1) - p;
F(2)=x(2)-tau;
