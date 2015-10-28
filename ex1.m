close all;
clear;

F=1000;
lambda = 0.2;
p=1e-3;
tp=1*lambda*sqrt(-2*log(p));
fp = (1/(pi*lambda))*sqrt((-log(p))/(2));

delta_t = 20*tp;

t=0:1/F:delta_t;

s_t = (1/(lambda*sqrt(2*pi)))*exp(-1*((t-tp).^2/(2*lambda^2)));

figure(1);
subplot(2,1,1);
hold on;
plot(t, s_t);
xlabel('Temps (t)');
ylabel('Signal');
title('Signaux s(t) et s_r(t)');
alpha = 0.8*exp(i*pi/4);
Tsec = 1;

sr_t = alpha*(1/(lambda*sqrt(2*pi)))*exp(-1*((t-tp-Tsec).^2/(2*lambda^2)));
plot(t,sr_t, 'r');

u_t = s_t + sr_t;
legend('s(t)','s_r(t)','Location','northeast')


subplot(2,1,2);
plot(t, u_t, 'g');
xlabel('Temps (t)');
ylabel('Signal');
title('Signal u(t)');


%%%%%%%%%%%%% FFT %%%%%%%%%%%%

fft_u = fft(u_t);
f=-F/2:1/delta_t:F/2;

figure(2);
subplot(2,1,1);

plot(f,abs(fft_u));
subplot(2,1,2);
plot(f, angle(fft_u));

figure(3);
fft_shift = fftshift(fft_u);

milieu = length(fft_shift)/2;

delta = (fp/F)*length(fft_shift);

subplot(2,1,1);
plot(f(milieu-delta:milieu+delta),abs(fft_shift(milieu-delta:milieu+delta)));
subplot(2,1,2);
plot(f(milieu-delta:milieu+delta),angle(fft_shift(milieu-delta:milieu+delta)));

%%%%%%%% FILTRE %%%%%%%%%%

w=1;
beta = 0.8;
val = (1-beta)/(2*w);

f_filtre = -fp:fp;
h_f = sqrt(w).*(abs(f_filtre)<=(val)) + sqrt((w/2)*(1+sin(pi*w/beta*((1/(2*w))-abs(f_filtre))))).*(abs(f_filtre) >= val & abs(f_filtre) <= val);

figure(4);
plot(f_filtre, abs(h_f));

figure(5);
plot(ifft(h_f));
