close all
for i = [6,10,20,40,100,200]
  a = load('acc.txt');
  v = load('vel.txt');
  d = load('dip.txt');
  padd = 2^16;
  a = [zeros(1,padd) a zeros(1,padd)];
  v = [zeros(1,padd) v zeros(1,padd)];
  d = [zeros(1,padd) d zeros(1,padd)];
  dt = 0.1;
  Fs = 1/dt;
  freq = 0:Fs/length(d):Fs/2;
  freq = freq(1:end-1);
  freq = freq*2*pi/0.05;
  FFTa = abs(fft(a)).^2;
  FFTa = FFTa(1:end/2);
  FFTa = FFTa/max(FFTa(freq<2));
  FFTd = abs(fft(d)).^2;
  FFTd = FFTd(1:end/2);
  FFTd = FFTd.*freq.^4;
  FFTd = FFTd/max(FFTd(freq<2));
  FFTv = abs(fft(v)).^2;
  FFTv = FFTv(1:end/2);
  FFTv = FFTv.*freq.^2;
  FFTv = FFTv/max(FFTv(freq<2));

  figure
  subplot(3,3,1)
  semilogy(freq,FFTd)
  title('EF Length Acceleration')
  axis([0,i,1e-8,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,2)
  semilogy(freq,FFTv)
  title('EF Velocity Acceleration')
  axis([0,i,1e-8,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,3)
  semilogy(freq,FFTa)
  title('EF Accleration Acceleration');
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  axis([0,i,1e-8,max(FFTa)])


  a = load('acc.txt');
  v = load('vel.txt');
  d = load('dip.txt');

  a = [zeros(1,padd) a zeros(1,padd)];
  v = [zeros(1,padd) v zeros(1,padd)];
  d = [zeros(1,padd) d zeros(1,padd)];
  dt = 0.1;
  Fs = 1/dt;
  freq = 0:Fs/length(d):Fs/2;
  freq = freq(1:end-1);
  freq = freq*2*pi/0.05;
  FFTa = abs(fft(a)).^2;
  FFTa = FFTa(1:end/2)./freq.^2;
  FFTa = FFTa/max(FFTa(and(freq>0.5,freq<2)));
  FFTd = abs(fft(d)).^2;
  FFTd = FFTd(1:end/2);
  FFTd = FFTd.*freq.^2;
  FFTd = FFTd/max(FFTd(freq<2));
  FFTv = abs(fft(v)).^2;
  FFTv = FFTv(1:end/2);
  FFTv = FFTv/max(FFTv(freq<2));

  subplot(3,3,4)
  semilogy(freq,FFTd)
  title('EF Length Velocity')
  axis([0,i,1e-10,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,5)
  semilogy(freq,FFTv)
  title('EF Velocity Velocity')
  axis([0,i,1e-10,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,6)
  semilogy(freq,FFTa)
  title('EF Acceleration Velocity');
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  axis([0,i,1e-10,max(FFTa)])

  a = load('acc.txt');
  v = load('vel.txt');
  d = load('dip.txt');

  a = [zeros(1,padd) a zeros(1,padd)];
  v = [zeros(1,padd) v zeros(1,padd)];
  d = [zeros(1,padd) d zeros(1,padd)];
  dt = 0.1;
  Fs = 1/dt;
  freq = 0:Fs/length(d):Fs/2;
  freq = freq(1:end-1);
  freq = freq*2*pi/0.05;
  FFTa = abs(fft(a)).^2;
  FFTa = FFTa(1:end/2)./freq.^4;
  FFTa = FFTa/max(FFTa(and(freq>0.5,freq<2)));
  FFTd = abs(fft(d)).^2;
  FFTd = FFTd(1:end/2);
  FFTd = FFTd/max(FFTd(freq<2));
  FFTv = abs(fft(v)).^2;
  FFTv = FFTv(1:end/2);
  FFTv = FFTv./freq.^2;
  FFTv = FFTv/max(FFTv(and(freq>0.5,freq<2)));

  subplot(3,3,7)
  semilogy(freq,FFTd)
  title('EF Length Dipole')
  axis([0,i,1e-14,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,8)
  semilogy(freq,FFTv)
  title('EF Velocity Dipole')
  axis([0,i,1e-14,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
  subplot(3,3,9)
  semilogy(freq,FFTa)
  title('EF Accleration Dipole');
  axis([0,i,1e-14,max(FFTa)])
  grid on
  hold on
  xline(0.5/0.05)
  hold off
 
end  