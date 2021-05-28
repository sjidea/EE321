%Sim_ASK_passband_noncoherent.m
clear, clf
b=1; M=2^b;
ss=[0;1];
Tb=1; Ts=b*Tb; Ns=40;
T=Ts/Ns; LB=4*Ns;
Es=2; A=sqrt(Es);
thrshld=A/2;

%passband ASK Waveform
wc=10*pi/Ts; wcT=wc*T;  t=[0:Ns-1]*T; wct=wc*t; tt=[0:LB-1]*T;
su=sqrt(2/Ts)*[cos(wct); sin(wct)]; suT=su*T;
sw(1,:)=zeros(1,Ns); sw(2,:)=A*su(1,:);


%% Noncoherent Envelope detector

SNRdBs=[1:10];  MaxIter=10000;           % Range of SNRdB and number of iterations
SNRbdBt=0:0.1:10;  SNRbt=10.^(SNRbdBt/10);

% make your own codes
load('message')
load('noise_i')
load('noise_q')
Ln = length(message);
LnSnr = length(SNRdBs);
my_t = [0:Ns*Ln-1]*T;

% message without noise
for i = 1:Ln 
    if message(i,1) == 1
        message(i,1) = ss(1,1);
    else
        message(i,1) = ss(2,1);
    end
end
r = zeros (1, Ns*Ln);          % message
for i = 1:Ns
    for k = 1:Ln
        r(1,i+(k-1)*Ns) = sw(message(k,1)+1,i);
    end
end
rm = r(1,Ns*Ln-LB+1:Ns*Ln);

%message with noise
n_i = zeros (1, Ns*Ln, LnSnr);
n_q = zeros (1, Ns*Ln, LnSnr);
n   = zeros (LnSnr, Ns*Ln);     % noise
rn  = zeros (LnSnr, Ns*Ln);     % message + noise
for i = 1:LnSnr
    for j = 1:Ln
        for k = 1:Ns
            n_i(1,k+Ns*(j-1),i) = noise_i(i,j,k);
            n_q(1,k+Ns*(j-1),i) = noise_q(i,j,k);
        end
    end
    n (i,1:Ns*Ln) = n_i(1,:,i).*cos(wc*my_t) - n_q(1,:,i).*sin(wc*my_t);
    rn(i,1:Ns*Ln) = r + n(i,1:Ns*Ln)*sqrt(Ns*(10^(-i/10)));
end
rnm = rn(10,Ns*Ln-LB+1:Ns*Ln);

%envelope output
y_noc = zeros(1,Ns*Ln);         % message * cos
y_nos = zeros(1,Ns*Ln);         % message * sin
y_Noc = zeros(LnSnr,Ns*Ln);     % (message + noise) * cos
y_Nos = zeros(LnSnr,Ns*Ln);     % (message + noise) * sin
for i = 1:Ln-1
    for h = 1:Ns
        ti = (i-1)*Ns+h;
        y_noc(1,ti) = sum( r(1,ti:ti+39).*suT(1,:) );
        y_nos(1,ti) = sum( r(1,ti:ti+39).*suT(2,:) );
        for j =1:10
            y_Noc(j,ti) = sum( rn(j,ti:ti+39).*suT(1,:) );
            y_Nos(j,ti) = sum( rn(j,ti:ti+39).*suT(2,:) );
        end
    end
end
y_no = y_noc.^2 + y_nos.^2;             % yc(t)^2 + ys(t)^2
y_No = y_Noc.^2 + y_Nos.^2;             % with noise

ymc  = zeros (1, LB);
yms  = zeros (1, LB);
ynmc = zeros (1, LB);
ynms = zeros (1, LB);
for i = 1:4
    for h = 1:Ns
        ti = (Ln-6+i)*Ns +h;
        ymc (1 ,(i-1)*Ns+h)    = sum(  r(1 ,ti+1:ti+40).*suT(1,:));
        yms (1 ,(i-1)*Ns+h)    = sum(  r(1 ,ti+1:ti+40).*suT(2,:));
        ynmc(10, (i-1)*Ns +h)  = sum( rn(10,ti+1:ti+40).*suT(1,:));
        ynms(10, (i-1)*Ns +h)  = sum( rn(10,ti+1:ti+40).*suT(2,:));
    end
end
ym  =  ymc.^2 +  yms.^2;            
ynm = ynmc.^2 + ynms.^2;

%envelope output samples at SNR=10
aa = [1 : Ln];
y_p = zeros(LnSnr,Ln);              % z^2(k*Ts)
for j = 1:LnSnr
    for i = 1:Ln
        y_p(j,i) = y_No(j,i*Ns-39);
    end
end
ab = [1:1000];
y_pm = y_p(10,1:1000);

%BER
y_pd  = zeros (LnSnr, Ln);          % D[k]
n_Ber = zeros (1,LnSnr);            % BER
for j = 1:LnSnr
    num = 0;
    for i = 1:Ln
        if sqrt(y_p(j,i)) > thrshld 
            y_pd(j,i) = 1;
        end
        if message(i) ~= y_pd(j,i)
            num = num +1;
        end
    end
    n_Ber(1,j) = num / Ln;
end
n_Bert = exp(-SNRbt/4)/2; 


subplot(221), plot(tt,rm,'black', tt, rnm, 'b:')
title('Recieived signal r(t)')
subplot(222), plot(tt, ym, 'black', tt, ynm , 'b:');
title('Envelope detector output')
subplot(223), plot(ab, y_pm, '*'); 
title('Envelope detector output samples at SNRdB=10dB')
subplot(224), semilogy(SNRbdBt, n_Bert, 'black', SNRdBs, n_Ber, 'b*');
title('Bit Error Rate for passband ASK Signaling')