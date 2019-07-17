clc;
clear all;
close all;

x = randi([0 1],1,10000);
tx = 2*x-1;
X12 = reshape(tx,2,5000);

snr = 0:0.2:25;

for i = 1:length(snr)
    error(i) = 0;
    h11 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    h12 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    h21 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    h22 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    
    n = 10^(-snr(i)/20).*randn(1,10000);
    n12 = reshape(n,2,5000);
    
    for k = 1:5000  
        H11 = [h11(k) h12(k); h21(k) h22(k)];
        Y(1:2,k) = H11*[X12(1,k); X12(2,k)] + [n12(1,k);n12(2,k)];  
    end   
    
    for m=1:5000
        H22 = [h11(m) h12(m); h21(m) h22(m)];  
        W = inv((ctranspose(H22)*H22))*ctranspose(H22);
        X_cap(1:2,m) = W*[Y(1,m); Y(2,m)];
    end
    
    X_capr = reshape(X_cap,1,10000);
    dec = real(X_capr)>0;
    
    for l = 1:10000
        if dec(l) ~= x(l)
            error(i) = error(i) +1;
        end    
    end
    
    p(i) = (1/2)*(1-sqrt(10^(snr(i)/10)/(10^(snr(i)/10)+1)));
    
end


semilogy(snr,(error/10000));
hold on;
semilogy(snr,p,'--','LineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR');
legend('Simulation','Theory')