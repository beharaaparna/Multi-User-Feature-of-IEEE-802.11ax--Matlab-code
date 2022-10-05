function F = root2d_new_fpa_ng_23_nonsat_dyn(x)
global W_ng;
% global Nra_ng;
global m_ul;
global m_limit;
global g;
global timeq;
global lambda;
global new_stas;
global dyn_Nra;
% global p;
Nra_ng=dyn_Nra;
Nra=Nra_ng;
sum=0;sum12=0;nra=new_stas;
Nra=Nra_ng;w=W_ng;r=Nra;mbar=m_ul;m=m_limit;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1_num=w*(1-((2*x(1))^(mbar+1)));
a1_den=1-(2*x(1));a1_den=2*r*a1_den;
a1=a1_num/a1_den;

a2_num=w*(2^mbar)/(2*r);

sumq=0;
for ii=0:1:mbar
    sumq=sumq+((x(1)/2)^ii);
end
a2=(1/w)*sumq;

sumq=0;
for ii=0:1:mbar
    sumq=sumq+((2*x(1))^ii);
end
a4=(w/(2*r))*sumq;

sumq=0;
for ii=mbar+1:1:m
    sumq=sumq+((x(1)^ii));
end
vv1=1/((2^mbar)*w);
vv2=(2^mbar*w)/(2*r);
a21=(vv1+vv2)*sumq;


% a3_num=((r)-2)/(2*r);
a3_num=((3*r)-1)/(2*r);
sumq1_num=1-(x(1)^(m+1));
sumq1_den=1-x(1);
sumq1=sumq1_num/sumq1_den;
a3=a3_num*sumq1;

sum1=0;
for i=mbar+1:1:m
    sum1=sum1+(x(1)^i);
end
a6=sum1/((2^mbar)*w);

ll=x(1)/2;
a7_num=1-(ll^(mbar+1));
a7_den=1-(ll);
a7=a7_num/a7_den;
a7=a7/w;


tau_den=a2+a3+a21+a4;
% tau_den=a1+a2+a3+a6+a7;
tau_num1=1-(x(1)^(m+1));
tau_num2=1-x(1);

den_val1=tau_den;
[delay1,delay2]=find_moments_mean_alt_drop(w,m_limit,mbar,x(1),Nra_ng);
% delay2=((m+1)*timeq)/((x(1)^(m+1)))
delay=(delay1+delay2)/2;
qq=min(1,lambda*delay*(10^(-6)));
rr=1-exp(-lambda*timeq*(10^(-6)));
den_val2=(1-qq)/rr;
if(x(1)~=0)
    den_val=den_val1+den_val2;
else
    qq=min(1,lambda*timeq*(10^(-6)));
    rr=1-exp(-lambda*timeq*(10^(-6)));
    den_val2=(1-qq)/rr;
    den_val=den_val1+den_val2;
end
tau_num=tau_num1/tau_num2;
tau=tau_num/den_val;

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
