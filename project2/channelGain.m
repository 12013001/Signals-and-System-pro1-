function [result] = channelGain(pilot)
dt = 1e-9; 
wc = 1e8;
R = 1000; % # of samples in 1us
N = 32;
fct = 1/dt;

%% ifft
ofdm = ifft(pilot);


%% add CP
cp = 4; % length of CP = 4
ofdm_cp = [ofdm((length(ofdm)-(cp-1)):(length(ofdm))),ofdm]; 

%% DAC and zero order hold
xz=reshape(repmat(ofdm_cp,R,1),1,[]);% dac and zero order hold
t = (1:length(xz))*dt;

%% modulation
xm = real(xz).*cos(2*pi*wc*t) + 1i*imag(xz).*sin(2*pi*wc*t);

%% channel
h = [0.5,zeros(1,1499),0.4,zeros(1,999),0.35,zeros(1,499), 0.3];
xc = filter(h,1,xm);

%% demodulation
xd = 2*(real(xc).*cos(2*pi*wc*t) + 1i*imag(xc).*sin(2*pi*wc*t));
[b,a] = butter(4,wc/(fct/2),'low');
xd_lpf = filter(b,a,xd); % LPF

%% ADC
xdd = reshape(xd_lpf,R,[]);
y_cp = mean(xdd,1);

%% remove cp
y = y_cp(cp+1:cp+32);
result = fft(y)./pilot/N; % derive the channel gain
end

