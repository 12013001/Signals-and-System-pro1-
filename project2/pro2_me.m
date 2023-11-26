dt = 1e-9;
fct = 1/dt; % CT sampling frequency
wc = 1e8; % carrier frequence
R = 1000; % # of samples in 1us
N = 32;
n = 0:31;
% generate signal
pilot = randn(1,32) + 1i*randn(1,32);
data = randn(1,32) + 1i*randn(1,32);
figure(1)
subplot(2,1,1),stem(n,real(data)),xlabel('n');title('real data to transmit')
subplot(2,1,2),stem(n,imag(data)),xlabel('n');title('imag data to transmit')

% ifft
ofdm = ifft(data);
figure(2)
subplot(2,1,1),stem(n,real(ofdm)),xlabel('n');title('real OFDM symbols')
subplot(2,1,2),stem(n,imag(ofdm)),xlabel('n');title('imag OFDM symbols')

% add CP
cp = 2; % length of CP
ofdm_cp = [ofdm((length(ofdm)-(cp-1)):(length(ofdm))),ofdm]; 
figure(3)
subplot(2,1,1),stem(0:length(ofdm_cp)-1,real(ofdm_cp)),xlabel('n');title('real OFDM symbols with CP')
subplot(2,1,2),stem(0:length(ofdm_cp)-1,imag(ofdm_cp)),xlabel('n');title('imag OFDM symbols with CP')

% DAC Transmitter
xa=upsample(ofdm_cp,R);% upsample the original signal
xz=reshape(repmat(ofdm_cp,R,1),1,[]);% zero order hold
t = (1:length(xz))*dt; % time vector
figure(4)
subplot(3,1,1),plot(t,real(xa)),xlabel('t'),title('pulse')
subplot(3,1,2),plot(t,real(xz)),xlabel('t'),title('real pulse shaping')
subplot(3,1,3),plot(t,imag(xz)),xlabel('t'),title('imag pulse shaping')

% RF-front end
xm = real(xz).*cos(2*pi*wc*t) + 1i*imag(xz).*sin(2*pi*wc*t);
x_fr = fft(xz)/length(xz);
x_frp = abs(fftshift(x_fr));
center = floor(length(x_frp)/2);
figure(5)
subplot(2,1,1),plot(t,real(xm)),ylabel('OFDM pilot RF'),xlabel('t');
subplot(2,1,2),stem(((-100):(100))./dt./length(x_fr)+1e8, x_frp((center-100):(center+100))); 
ylabel('freq of OFDM pilot RF');

% channel
h = [0.5,zeros(1,1499),0.4,zeros(1,999),0.35,zeros(1,499), 0.3];
xc = filter(h,1,xm);
figure(6)
subplot(2,1,1), plot((1:length(xc)).*dt, real(xm)), title('real part before channel');
subplot(2,1,2), plot((1:length(xc)).*dt, real(xc)), title('real part after channel');

% demodulation
xd = 2*(real(xc).*cos(2*pi*wc*t) + 1i*imag(xc).*sin(2*pi*wc*t));
[b,a] = butter(4,wc/(fct/2),'low');
xd_lpf = filter(b,a,xd); % LPF
figure(7)
subplot(2,1,1),plot(t,real(xd_lpf)),xlabel('t'),title('real part of the demodulated signal')
subplot(2,1,2),plot(t,imag(xd_lpf)),xlabel('t'),title('imag part of the demodulated signal')

% ADC
xdd = reshape(xd_lpf,R,[]);
y_cp = mean(xdd,1);
figure(8)
subplot(2,1,1),stem((1:length(y_cp)),real(y_cp)),title('real part'); 
subplot(2,1,2),stem((1:length(y_cp)),imag(y_cp)), title('imag part');

% remove cp
y = y_cp(cp+1:cp+32);
figure(9)
subplot(2,1,1),stem(real(y)),title('real part'); 
subplot(2,1,2),stem(imag(y)), title('imag part');

%  remove the channel effect
Xc = fft(y)./channelGain(pilot)/N;
xc = ifft(Xc);

ofdm_fft = fft(y);
figure(10)
subplot(4,1,1),stem(n,real(ofdm_fft)),xlabel('n');title('real OFDM symbols received')
subplot(4,1,2),stem(n,real(ofdm)),xlabel('n');title('real OFDM symbols')
subplot(4,1,3),stem(n,imag(ofdm_fft)),xlabel('n');title('imag OFDM symbols received')
subplot(4,1,4),stem(n,imag(ofdm)),xlabel('n');title('imag OFDM symbols')

figure(11)
subplot(4,1,1), stem(0:31, real(ofdm_fft)),title(' real OFDM received' ),xlabel(' t' )
subplot(4,1,2), stem(0:31,real (data)),title(' real OFDM pilot'), xlabel('t')
subplot(4,1,3), stem(0:31, imag(ofdm_fft)), title(' imag OFDM received' ),xlabel(' t' )
subplot(4,1,4), stem(0:31, imag(data)), title(' imag OFDM pilot'), xlabel('t')


channelGain=channelGain(pilot);
channelGain_timedomain=ifft(channelGain);
figure(12)
subplot(2,1,1),stem(0:length(channelGain)-1,real(channelGain)),title('channelGain in frequency domain')
subplot(2,1,2),stem(0:length(channelGain)-1,real(channelGain_timedomain)),title('channelGain in time domain')

figure(13)
stem(0:length(xc)-1,real(xc)),xlabel('n'),title('received siganl')

% difference
figure(14)
plot(0:length(ofdm)-1,real(ofdm),'r-')
hold on
plot(0:length(y)-1,real(xc),'b--')
title('length of cp = 2')
legend('the original singal','the received signal')
hold off;
