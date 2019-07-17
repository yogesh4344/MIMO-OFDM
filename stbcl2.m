clc;
clear all;
close all;

N = 64*10000;
x = randi([0 1],1,N);
tx = 2*x-1;
xr = reshape(tx,64,N/64);
ift = ifft(xr);
X12 = (1/sqrt(2)).*reshape(ift,2,N/2);

snr = 0:1:15;

for i = 1:length(snr)
    error(i) = 0;
    error1(i) = 0;
    h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
    
    n1 = 10^(-snr(i)/20).*randn(4,N/2);
    
    for k = 1:N/2  
        H11 = [h11(k) h12(k); h21(k) h22(k); conj(h12(k)) -conj(h11(k)); conj(h22(k)) -conj(h21(k))];
        Y(1:4,k) = H11*[X12(1,k); X12(2,k)] + [n1(:,k)];  
    end   
    
    for m=1:N/2
        H22 = [h11(m) h12(m); h21(m) h22(m); conj(h12(m)) -conj(h11(m)); conj(h22(m)) -conj(h21(m))];
        W = inv((ctranspose(H22)*H22))*ctranspose(H22);
        X_cap(1:2,m) = W*[Y(1,m); Y(2,m); Y(3,m); Y(4,m)];
    end
    
    X_capr = fft(reshape(X_cap,64,N/64));
    X_capr1 = reshape(X_capr,1,N);
    dec = real(X_capr1)>0;
    
    for l = 1:N
        if dec(l) ~= x(l)
            error(i) = error(i) +1;
        end    
    end
end


semilogy(snr,(error/N),'black');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for 2X2 MIMO using Alamouti(STBC)');