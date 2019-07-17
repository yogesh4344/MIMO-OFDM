clc;
clear;
close all;
            %OFDM of 16 points FFT and IFFT

             %for 1X2 MIMO setup
             f = [16 64];
             
             for z = 1:length(f)
N = 256*100;
x = randi([0 1],1,N);
tx = 2*x-1;
x1 = reshape(tx,1,N/1);
x11 = reshape(x1(1,:),N/(1*f(z)),f(z));
ift1 = ifft(x11);
X12 = [reshape(ift1,1,N)];
            snr = 0:1:40;
            for i = 1:length(snr)
                error(i) = 0;
                h11 = (1/sqrt(2))*(randn(1,N) + j*randn(1,N));
                h21 = (1/sqrt(2))*(randn(1,N) + j*randn(1,N));
           

                n = 10^(-snr(i)/20).*randn(1,N);
                n12 = reshape(n,1,N);

                for k = 1:N 
                    H12 = [h11(k);h21(k)];
                    Y(1:2,k) = H12*X12(k) + n12(k);  
                end   

                for m=1:N
                    H12_m = [h11(m);h21(m)];  
                    W = inv((ctranspose(H12_m)*H12_m))*ctranspose(H12_m);
                    X_cap(1,m) = W*[Y(1,m); Y(2,m)];
                end

    x_c = fft(reshape(X_cap(1,:),N/(1*f(z)),f(z)));
    X_capr1 = [reshape(x_c,1,N)];
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
  
            
            
            
            %for 2X2 MIMO setup
N = 256*100;
x = randi([0 1],1,N);
tx = 2*x-1;
x1 = reshape(tx,2,N/2);
x11 = reshape(x1(1,:),N/(2*f(z)),f(z));
ift1 = ifft(x11);
x12 = reshape(x1(2,:),N/(2*f(z)),f(z));
ift2 = ifft(x12);
X22 = [reshape(ift1,1,N/2); reshape(ift2,1,N/2)];
            snr = 0:1:40;
              for i = 1:length(snr)
                error(i) = 0;
                h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));

                n = 10^(-snr(i)/20).*randn(1,N);
                n22 = reshape(n,2,N/2);

                for k = 1:N/2  
                    H22 = [h11(k) h12(k); h21(k) h22(k)];
                    Y(1:2,k) = H22*[X22(1,k); X22(2,k)] + [n22(1,k);n22(2,k)];  
                end   

                for m=1:N/2
                    H22_m = [h11(m) h12(m); h21(m) h22(m)];  
                    W = inv((ctranspose(H22_m)*H22_m))*ctranspose(H22_m);
                    X_cap1(1:2,m) = W*[Y(1,m); Y(2,m)];
                end

    x_c = fft(reshape(X_cap1(1,:),N/(2*f(z)),f(z)));
    x_c2 = fft(reshape(X_cap1(2,:),N/(2*f(z)),f(z)));
    X_capr1 = [reshape(x_c,1,N/2); reshape(x_c2,1,N/2)];
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
              
              
              
% % %             %for 2X3 MIMO setup
N = 256*100;
x = randi([0 1],1,N);
tx = 2*x-1;
x1 = reshape(tx,2,N/2);
x11 = reshape(x1(1,:),N/(2*f(z)),f(z));
ift1 = ifft(x11);
x12 = reshape(x1(2,:),N/(2*f(z)),f(z));
ift2 = ifft(x12);
X23 = [reshape(ift1,1,N/2); reshape(ift2,1,N/2)];
            snr = 0:1:40;
              for i = 1:length(snr)
                error(i) = 0;
                h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h31=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h32=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                
                
                n = 10^(-snr(i)/20).*randn(1,N);
                n23 = 10^(-snr(i)/20).*randn(3,N/2);

                for k = 1:N/2  
                    H23 = [h11(k) h12(k); h21(k) h22(k); h31(k) h32(k)];
                    Y(1:3,k) = H23*[X23(1,k); X23(2,k)] + [n23(1,k);n23(2,k);n23(3,k)];  
                end   

                for m=1:N/2
                    H23_m = [h11(m) h12(m); h21(m) h22(m);h31(m) h32(m)];  
                    W = inv((ctranspose(H23_m)*H23_m))*ctranspose(H23_m);
                    X_cap3(1:2,m) = W*[Y(1,m); Y(2,m); Y(3,m)];
                end

    x_c = fft(reshape(X_cap3(1,:),N/(2*f(z)),f(z)));
    x_c2 = fft(reshape(X_cap3(2,:),N/(2*f(z)),f(z)));
    X_capr1 = [reshape(x_c,1,N/2); reshape(x_c2,1,N/2)];
                X_capr = reshape(X_capr1,1,N);
                dec = real(X_capr)>0;

                for l = 1:N
                    if dec(l) ~= x(l)
                        error(i) = error(i) +1;
                    end    
                end

              end

semilogy(snr,(error/N));
              
              
              
% % %             %for 2X4 MIMO setup
N = 256*100;
x = randi([0 1],1,N);
tx = 2*x-1;
x1 = reshape(tx,2,N/2);
x11 = reshape(x1(1,:),N/(2*f(z)),f(z));
ift1 = ifft(x11);
x12 = reshape(x1(2,:),N/(2*f(z)),f(z));
ift2 = ifft(x12);
X23 = [reshape(ift1,1,N/2); reshape(ift2,1,N/2)];
            snr = 0:1:40;
              for i = 1:length(snr)
                error(i) = 0;
                h11 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h12 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h21 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h22 = (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h31=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h32=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h41=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                h42=  (1/sqrt(2))*(randn(1,N/2) + j*randn(1,N/2));
                
                
                
                n = 10^(-snr(i)/20).*randn(1,2*N);
                n24 = reshape(n,4,N/2);

                for k = 1:N/2  
                    H24 = [h11(k) h12(k); h21(k) h22(k); h31(k) h32(k); h41(k) h42(k)];
                    Y(1:4,k) = H24*[X23(1,k); X23(2,k)] + [n24(1,k);n24(2,k);n24(3,k);n24(4,k)];  
                end   

                for m=1:N/2
                    H24_m = [h11(m) h12(m); h21(m) h22(m);h31(m) h32(m);h41(m) h42(m)];  
                    W = inv((ctranspose(H24_m)*H24_m))*ctranspose(H24_m);
                    X_cap4(1:2,m) = W*[Y(1,m); Y(2,m); Y(3,m);Y(4,m)];
                end

    x_c = fft(reshape(X_cap4(1,:),N/(2*f(z)),f(z)));
    x_c2 = fft(reshape(X_cap4(2,:),N/(2*f(z)),f(z)));
    X_capr1 = [reshape(x_c,1,N/2); reshape(x_c2,1,N/2)];
                X_capr = reshape(X_capr1,1,N);
                dec = real(X_capr)>0;

                for l = 1:N
                    if dec(l) ~= x(l)
                        error(i) = error(i) +1;
                    end    
                end

              end
            semilogy(snr,(error/N));
             end
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR');
%legend('1X2','2X2','2x3','2X4');
legend('16 2X2','16 2x3','16 2X4','64 2X2','64 2x3','64 2X4');