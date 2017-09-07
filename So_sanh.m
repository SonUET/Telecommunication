clear all
N = 10^6;
Eb_N0_dB = [0:30];

% BPSK Modulation
ipBPSK = rand(1,N)>0.5; % tao cac bit 0, 1 voi xac suat bang nhau
s1 = 2*ipBPSK-1; % Thuc hien dieu che: 0 -> -1; 1 -> 1 
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % chuan hoa vecto Gauss co nang luong trung binh = 1 hay 0dB 
for ii = 1:length(Eb_N0_dB)
   % Nhieu
   y1 = s1 + 10^(-Eb_N0_dB(ii)/20)*n; % AWGN

   % Bo thu - hard decision decoding
   ipHat1 = real(y1)>0;

   % Dem so loi
   nErr1(ii) = size(find([ipBPSK- ipHat1]),2);
end
simBer1 = nErr1/N; % BER theo mo phong
theoryBer1 = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % BER theo ly thuyet


% QPSK Modulation
ipHat2 = zeros(1,N);
for ii = 1:length(Eb_N0_dB)
ipQPSK = (2*(rand(1,N)>0.5)-1) + j*(2*(rand(1,N)>0.5)-1); 
s2 = (1/sqrt(2))*ipQPSK; % normalization of energy to 1
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

y2 = s2 + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

% demodulation
y_re = real(y2); % real
y_im = imag(y2); % imaginary
ipHat2(find(y_re < 0 & y_im < 0)) = -1 + -1*j;
ipHat2(find(y_re >= 0 & y_im > 0)) = 1 + 1*j;
ipHat2(find(y_re < 0 & y_im >= 0)) = -1 + 1*j;
ipHat2(find(y_re >= 0 & y_im < 0)) = 1 - 1*j;

nErr2(ii) = size(find([ipQPSK- ipHat2]),2); % couting the number of errors
end
simSer_QPSK = nErr2/N;
theorySer_QPSK = erfc(sqrt(0.5*(10.^(Eb_N0_dB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(Eb_N0_dB/10))))).^2;

% 16-PSK Modulation
M = 16;
thetaMpsk = [0:M-1]*2*pi/M;
for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    ipPhase = randsrc(1,N,thetaMpsk);
    ip3 = exp(j*ipPhase);   
    s3 = ip3; % normalization of energy to 1
    
    % noise
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance 
    
    y3 = s3 + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise 
     
    % demodulation
    % finding the phase from [-pi to +pi]
    opPhase = angle(y3); 
    % unwrapping the phase i.e. phase less than 0 are 
    % added 2pi
    opPhase(find(opPhase<0)) = opPhase(find(opPhase<0)) + 2*pi;

    % rounding the received phase to the closest 
    % constellation
    ipHat2 = 2*pi/M*round(opPhase/(2*pi/M))	;
    % as there is phase ambiguity for phase = 0 and 2*pi,
    % changing all phases reported as 2*pi to 0.
    % this is to enable comparison with the transmitted phase
    ipHat2(find(ipHat2==2*pi)) = 0; 

    % counting errors
    nErr3(ii) = size(find([ipPhase- ipHat2]),2); % couting the number of errors
end
simBer3 = nErr3/N;
theoryBer3 = erfc(sqrt(10.^(Eb_N0_dB/10))*sin(pi/M));

%16-QAM Modulation
alpha16qam = [-3 -1 1 3]; % 16-QAM alphabets
for ii = 1:length(Eb_N0_dB)
    ipQAM = randsrc(1,N,alpha16qam) + j*randsrc(1,N,alpha16qam);
    s4 = (1/sqrt(10))*ipQAM; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

    y4 = s4 + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    y_re2 = real(y4); % real part
    y_im2 = imag(y4); % imaginary part

    ipHat_re(find(y_re2< -2/sqrt(10)))           = -3;
    ipHat_re(find(y_re2 > 2/sqrt(10)))           =  3;
    ipHat_re(find(y_re2>-2/sqrt(10) & y_re2<=0))  = -1;
    ipHat_re(find(y_re2>0 & y_re2<=2/sqrt(10)))   =  1;

    ipHat_im(find(y_im2< -2/sqrt(10)))           = -3;
    ipHat_im(find(y_im2 > 2/sqrt(10)))           =  3;
    ipHat_im(find(y_im2>-2/sqrt(10) & y_im2<=0))  = -1;
    ipHat_im(find(y_im2>0 & y_im2<=2/sqrt(10)))   =  1;
    ipHat2 = ipHat_re + j*ipHat_im;
    nErr4(ii) = size(find([ipQAM- ipHat2]),2); % couting the number of errors
end

simBer4 = nErr4/N;
theoryBer4 = 3/2*erfc(sqrt(0.1*(10.^(Eb_N0_dB/10))));

%64-QAM Modulation
M2 = 64;
k = sqrt(1/((2/3)*(M2-1))); % normalizing factor
m = [1:sqrt(M2)/2]; % alphabets
alphaMqam = [-(2*m-1) 2*m-1];
for ii = 1:length(Eb_N0_dB)

    ip5 = randsrc(1,N,alphaMqam) + j*randsrc(1,N,alphaMqam);
    s5 = k*ip5; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

    y5 = s5 + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    y_re3 = real(y5)/k; % real part
    y_im3 = imag(y5)/k; % imaginary part
    
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    ipHat_re2 = 2*floor(y_re3/2)+1;
    ipHat_re2(find(ipHat_re2>max(alphaMqam))) = max(alphaMqam);
    ipHat_re2(find(ipHat_re2<min(alphaMqam))) = min(alphaMqam);
	     
    % rounding to the nearest alphabet
    % 0 to 2 --> 1
    % 2 to 4 --> 3
    % 4 to 6 --> 5 etc
    ipHat_im2 = 2*floor(y_im3/2)+1;
    ipHat_im2(find(ipHat_im2>max(alphaMqam))) = max(alphaMqam);
    ipHat_im2(find(ipHat_im2<min(alphaMqam))) = min(alphaMqam);
    
    ipHat2 = ipHat_re2 + j*ipHat_im2; 
    nErr5(ii) = size(find([ip5- ipHat2]),2); % counting the number of errors
end
simSer = nErr5/N;
theorySer = 2*(1-1/sqrt(M2))*erfc(k*sqrt((10.^(Eb_N0_dB/10)))) - (1-2/sqrt(M2) + 1/M2)*(erfc(k*sqrt((10.^(Eb_N0_dB/10))))).^2;


% plot
close all
figure
semilogy(Eb_N0_dB,theoryBer1,'LineWidth', 2);
hold on
semilogy(Eb_N0_dB,simBer1,'*');
semilogy(Eb_N0_dB,theorySer_QPSK,'LineWidth', 2);
semilogy(Eb_N0_dB,simSer_QPSK,'*');
semilogy(Eb_N0_dB,theoryBer3,'LineWidth', 2);
semilogy(Eb_N0_dB,simBer3,'*');
semilogy(Eb_N0_dB,theoryBer4,'LineWidth', 2);
semilogy(Eb_N0_dB,simBer4,'*');
semilogy(Eb_N0_dB,theorySer,'LineWidth', 2);
semilogy(Eb_N0_dB,simSer,'*');
axis([0 30 10^-6 1])
grid on
legend('ly thuyet-BPSK', 'mo phong-BPSK','ly thuyet-QPSK', 'mo phong-QPSK', 'ly thuyet-16 PSK', 'mo phong-16 PSK',...
'ly thuyet-16 QAM', 'mo phong-16 QAM','ly thuyet-64 QAM','mo phong-64 QAM');
xlabel('Eb/No, dB');
ylabel('Symbol Error Rate');
title('Symbol error probability curve')
