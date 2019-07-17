clc;
clear all;
close all;
f=[16 32 64 128];
for z=1:4
N = 256*10;
x = randi([0 1],1,N);
tx = 2*x-1;
x1 = reshape(tx,2,N/2);
x11 = reshape(x1(1,:),N/(2*f(z)),f(z));
ift1 = ifft(x11);
x12 = reshape(x1(2,:),N/(2*f(z)),f(z));
ift2 = ifft(x12);
X12 = [reshape(ift1,1,N/2); reshape(ift2,1,N/2)];
% Xr = reshape(tx,f(z),N/f(z));
% ift = ifft(Xr);
% X12 = reshape(ift,2,N/2);

snr = 0:0.2:25;

for i = 1:length(snr)
    error(i) = 0;
    h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    
    n = 10^(-snr(i)/20).*randn(1,N);
    n12 = reshape(n,2,N/2);
    
    for k = 1:N/2  
        H11 = [h11(k) h12(k); h21(k) h22(k)];
        Y(1:2,k) = H11*[X12(1,k); X12(2,k)] + [n12(1,k);n12(2,k)];  
    end   
    
    for m=1:N/2
        H22 = [h11(m) h12(m); h21(m) h22(m)];  
        W = inv((ctranspose(H22)*H22))*ctranspose(H22);
        X_cap(1:2,m) = W*[Y(1,m); Y(2,m)];
    end
    
    
    x_c = fft(reshape(X_cap(1,:),N/(2*f(z)),f(z)));
    x_c2 = fft(reshape(X_cap(2,:),N/(2*f(z)),f(z)));
    X_capr1 = [reshape(x_c,1,N/2); reshape(x_c2,1,N/2)];
%     X_fft = fft(reshape(X_cap,f(z),N/f(z)));
    X_capr = reshape(X_capr1,1,N);
    dec = real(X_capr)>0;
    
    for l = 1:N
        if dec(l) ~= x(l)
            error(i) = error(i) +1;
        end    
    end
    
end


semilogy(snr,(error/N));
hold on;
end
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR');
legend(' 16 fft mimo-ofdm',' 32 fft mimo-ofdm','64 fft mimo-ofdm','128 fft mimo-ofdm');
