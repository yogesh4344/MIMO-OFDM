clc;
clear all;
close all;

x = randi([0 1],1,10000);
tx = 2*x-1;
X12 = reshape(tx,2,5000);

snr = 0:0.1:15;

for i = 1:length(snr)
    error(i) = 0;
    h1 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    h2 = (1/sqrt(2))*(randn(1,5000) + j*randn(1,5000));
    
    n1 = 10^(-snr(i)/20).*randn(1,5000);
    n2 = 10^(-snr(i)/20).*randn(1,5000);
    
    for k = 1:5000  
        H11 = [h1(k) h2(k); conj(h2(k)) -conj(h1(k))];
        Y(1:2,k) = H11*[X12(1,k); X12(2,k)] + [n1(k);n2(k)];  
    end   
    
    for m=1:5000
        H22 = [h1(m) h2(m); conj(h2(m)) -conj(h1(m))]; 
        %W = inv((ctranspose(H22)*H22))*ctranspose(H22);
        X_cap(1:2,m) = inv(H22)*[Y(1,m); Y(2,m)];
    end
    
    X_capr = reshape(X_cap,1,10000);
    dec = real(X_capr)>0;
    
    for l = 1:10000
        if dec(l) ~= x(l)
            error(i) = error(i) +1;
        end    
    end

    p1 = 0.5 - 0.5*(1+(2/10^(snr(i)/10)))^(-0.5);
    p(i) = (p1^2)*(1+2*(1-p1));
    pm = 0.5 - 0.5*(1+(1/10^(snr(i)/10)))^(-0.5);
    pm1(i) = (pm^2)*(1+2*(1-pm));
end


semilogy(snr,2.*(error/10000));
hold on;
semilogy(snr,p,'--','LineWidth',2);
semilogy(snr,pm1,'--','LineWidth',2);
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR for 2x1 MIMO');
legend('Simulation','Theory','mrc')