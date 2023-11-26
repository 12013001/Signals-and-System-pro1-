clear;
%DAC
input=[1];
inX_p=[1,zeros(1,999)];
inX_c=input*ones(1,1000);
n=1:length(inX_p);
t=n/1000;

%Actual wireless channel h(t)
h=[0.5, 0, 0, 0.4, 0, 0.35, 0.3, 0];
ht=upsample(h,1000/2);%dt=1/2*T
t_ht=(1:length(ht))/1000;
CR=conv(inX_c,ht);%Channel receive

%ADC
CR=[CR,0];
n1=1:length(CR);
t1=n1/1000;
A=reshape(CR,1000,[]);
hn=mean(A,1);

figure(1);
stem(-4:4,[0,0,0,0,input,0,0,0,0],'filled');
title('Impluse');
figure(2);
stem(t,inX_p);
title('x_p');xlim([-2,4]);
figure(3);
plot(t,inX_c);
title('x_c');xlim([0,4]);
figure(4);
plot(t_ht,ht);
title('h(t)');
figure(5);
plot(t1,CR);
title('Channel received');
figure(6);
stem(0:length(hn)-1,hn);
title('h[n]');