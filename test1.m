clc;
clear all;
close all;
N = 16*10^3;
x = randi([0 1],1,N);
tx = 2*x-1;
Xr = reshape(tx,16,N/16);
ift = ifft(Xr);
X22 = reshape(ift,2,N/2);

snr = 0:0.2:25;

for i = 1:length(snr)
    error(i) = 0;
    h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h13 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h14 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h23 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h24 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    
    n = 10^(-snr(i)/20).*randn(1,N);
    n22 = reshape(n,2,N/2);
    
    for k = 1:N/2  
        H22 = [h11(k) h12(k); h21(k) h22(k)];
        Y(1:2,k) = H22*[X22(1,k); X22(2,k)] + [n22(1,k);n22(2,k)];  
    end   
    
    for m=1:N/2
        D22 = [h11(m) h12(m); h21(m) h22(m)];  
        W = inv((ctranspose(D22)*D22))*ctranspose(D22);
        X_cap(1:2,m) = W*[Y(1,m); Y(2,m)];
    end
    
    X_fft = fft(reshape(X_cap,16,N/16));
    X_capr = reshape(X_fft,1,N);
    dec = real(X_capr)>0;
    
    for l = 1:N
        if dec(l) ~= x(l)
            error(i) = error(i) +1;
        end    
    end
    
end


semilogy(snr,(error/N));
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR');
legend(' 16 fft mimo-ofdm')