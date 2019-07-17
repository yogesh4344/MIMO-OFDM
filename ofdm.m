clc;
clear all;
close all;


x = randi([0 1],1,64*100000);
tx = 2*x-1;
xr = reshape(x,64,100000);

snr = 0:0.2:25;

for i = 1:length(snr)
    error(i) = 0;
    ift = ifft(xr);
    y = awgn(ift,snr(i));
    ft = fft(y);
    rx = reshape(ft,1,64*100000);
    rx1 = real(rx)>0;
    
    for l = 1:100000
     if rx1(l) ~= x(l)
       error(i) = error(i) +1;
     end    
    end
    
end

semilogy(snr,(error/(64*100000)));
