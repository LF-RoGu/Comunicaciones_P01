%% P01_705694

close all; clc; clear;

filename= 'Never_Gonna.wav';

fir_order = 100;
Fs_audio = 44100;
Ts_audio = 1/Fs_audio;
t = 10; %t = 10 seconds

samples = [1,t*Fs_audio];

[x,Fs] = audioread(filename,samples);   %obtencion del archivo wav

info = audioinfo(filename);

info_bps = info.BitsPerSample; %info_bps (info_bits per sample)

    %Signal to Mono
xMono = (x(:,1)+ x(:,2))/2;

%% Filtrar la señal 
    %LPF Comunication Channel  
B = 15e3;
Fc = B / (Fs_audio/2); %Normalized Freq 0.6803

fc_15KHz = [0 Fc Fc 1];
m = [1 1 0 0];
order = fir_order;
LPF = fir2(order, fc_15KHz, m); % LPF (Low Pass Filter)
%fvtool(LPF, 'Analysis', 'impulse');
xa = filter(LPF, 1, xMono); %señal analogica filtrada
    % If mean it is equal to ZERO
if 1 == mean(xa)
    % Then we can assume that the variation is the P
    var(xa)
end
P_signal = var(xa); %verificar la potencia 
%% Normalize the Power
    % Calculate the power
pot = sum(xa.*xa)/numel(xa);
xa = xa/sqrt(pot);
    % check if the power of the signal
    % If mean it is equal to ZERO
if 0 == mean(xa)
   % Then we can assume that the variation is the P
   var(xa)
end
    % At this point must be equal to 1
P_signal = var(xa); %verificar la potencia 
%% SNR
N0 = 1./(B.*10.^(0:0.3:3)); % Vector PSD del ruido
P_noise = B*N0; % Vector de Pot del Ruido Filtrado
P_noise_dB = 10.*log10(P_noise); % Pot. Ruido en Decibeles
SNR_A = P_signal ./ P_noise; % Relacion Señal a Ruido
SNR_A_dB = 10*log10(SNR_A); % SNR en dB

%% AWGN (White noise all frequency)
for i = 1 : numel(N0)
    % noise samples
    noise = sqrt(P_noise(i)) * randn(1, numel(xa)); 
    
        % Add up noise to the signal
    xa_noised = xa + noise';
        % Normilized the signal between 1 & -1
    xa_noised = xa_noised ./ max(abs(xa_noised));

    filename = strcat('AnalogSignal_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
    audiowrite(filename, xa_noised, Fs); %obtencion del archivo 
end  
%% Reproduce sound
%soundsc(xMono,Fs_audio);

%% Digital Transmition
    %Source: https://www.mathworks.com/matlabcentral/answers/75379-how-to-convert-an-input-sine-wave-into-an-8-bit-digital-signal
    %Source: https://www.mathworks.com/matlabcentral/answers/103315-how-to-convert-digital-data-into-analog-data-using-matlab-code
    % Convert values from xMono to bin

swing = (2^info_bps-1)/2; 
xMono2bi = round(xMono*swing+swing);
    %For Default it is the 'right-msb'
b = de2bi(xMono2bi,info_bps,'left-msb');
    %Transponse operation
b = b';

bits_tx = b(:); %Concatena el resto de bits, para que sea un solo vector

%% 
    %This goes if beta it is equal to ZERO
%B = 15e3;
%Rb_max = B2*2;
%% 
beta = 0.35;
D = 6;
Fs = 96e3;
Rb = (2*B) / (1+beta);
Ts = 1/Fs;
mp = ceil(Fs/Rb); %mp = 5;
Tp = mp*Ts;
energy = Tp;

[p,t] = rcpulse(beta, D, Tp, Ts, 'srrc', energy);

p_energy = Ts*sum(p.*p);

p = p/sqrt(p_energy); 

%wvtool(p);

    %Polar Signal NRZ
PBP = p; % PSPB (Polar Base Pulse)

PS_s = bits_tx;

PS_s( PS_s == 0 ) = -1;

s = zeros(1,numel(PS_s)*mp);

s(1:mp:end) = PS_s;

PSLC = conv(PBP, s) / mp; %PSLC (Polar Signal Line Code)
    %Normilize pulse
PSLC = PSLC/sqrt( sum(PSLC.*PSLC)/numel(PSLC) );
P_PSLC_d = var(PSLC);


%% SNR Digital 
N0 = 1./(B.*10.^(0:0.3:3)); % Vector PSD del ruido
P_noise_d = B*N0; % Vector de Pot del Ruido Filtrado
P_noise_dB_d = 10.*log10(P_noise_d); % Pot. Ruido en Decibeles
SNR_D = P_PSLC_d ./ P_noise_d; % Relacion Señal a Ruido
SNR_D_dB = 10*log10(SNR_D); % SNR en dB

%% Add noise to Line Code 

% AWGN (White noise all frequency) 
%generate noise vector and add to polar NRZ line code with SRRC as pulse
%base

%match filter
%for the match filter we used the same form of the base pulse (SRRC)

[h_srrc t] = rcpulse(beta,D,Tp,Ts,'srrc',energy);
hsrrc_energy = Ts*sum(h_srrc.*h_srrc);
h_srrc = h_srrc/sqrt(hsrrc_energy); %normalized power of the filter

for i = 1 : numel(N0)
    % noise samples
    noise = sqrt(P_noise(i)) * randn(1, numel(PSLC)); %created noise vector
    
    % Add up noise to the signal
    PSLC_noised = PSLC + noise;
    
    % filtered line code + noise with the 15kh LPF

%figure; stem(PSLC(1:mp*10));

    PRS_rx = conv(LPF,PSLC_noised); %PSRS (Polar Received Signal)
%%figure; pwelch(PRS_rx,500,300,500,Fs,'power'); title('Polar Fc 0.6803');

%%figure; plot(PSLC(1:mp*20)); hold on; plot(PRS_rx(1:mp*20));

    %filtered the recieved line code in the match filter
    match_filter_recived = conv(h_srrc,PRS_rx);
    
    %make the sampling and desition of the line code recived
 % Tenemos order/2 + y - por el retardo del order del filtro
% Take a sample each mp/2 in time domain to recover the data
    rx_match_filter = match_filter_recived( ( mp/2 + ((order/2)+(numel(p)/2)) : mp : (end - ((order/2)+(numel(p)/2)))));

%en este punto, se obtiene la información recivada y pasada por el filtro
%que emula el canal y el match filter (recepcion) 
%% 
    audio_rx = sign(rx_match_filter); %Obtenemos valores entre 1 y -1 del arreglo

    bits_rx = (audio_rx+1)/2; %Normalizamos para que los valores sean 1 o 0

    bits_rx = bits_rx(4:end - 3);
    bits_rx = bits_rx';

    error = (sum( xor(bits_tx,bits_rx) ) / numel(bits_rx))*100;

%Convert vector to mat
    mat = vec2mat(bits_rx,16);
%Convert mat to dec
    audio_bin2dec = bi2de(mat,'left-msb');

    audio_mat = vec2mat(audio_bin2dec, 1 , Fs_audio);

%Normilize audio between 1 & -1
    %audio_mat = 2*mat2gray(audio_mat) - 1;
    audio_file = normalize(audio_mat, 'range', [-1 1]);

    filename = strcat('AnalogSignalRx_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
    audiowrite(filename, audio_file, Fs_audio); %obtencion del archivo 

end

%%
%figure; plot(PSLC_noised(1:mp*50));
%figure; plot(match_filter_recived(1:mp*50));
%figure; plot(rx_match_filter(1:mp*50));

%% Eye Pattern

EP = comm.EyeDiagram('SampleRate',Fs*mp,'SamplesPerSymbol',mp);

match_filter_recived = match_filter_recived';

%Eye Pattern Received Signal
EP(match_filter_recived);

%% 

audio_rx = sign(rx_match_filter); %Obtenemos valores entre 1 y -1 del arreglo

bits_rx = (audio_rx+1)/2; %Normalizamos para que los valores sean 1 o 0

bits_rx = bits_rx(4:end - 3);
bits_rx = bits_rx';

error = (sum( xor(bits_tx,bits_rx) ) / numel(bits_rx))*100;

%Convert vector to mat
mat = vec2mat(bits_rx,16);
%Convert mat to dec
audio_bin2dec = bi2de(mat,'left-msb');

audio_mat = vec2mat(audio_bin2dec, 1 , Fs_audio);

%Normilize audio between 1 & -1
%audio_mat = 2*mat2gray(audio_mat) - 1;
audio_mat = normalize(audio_mat, 'range', [-1 1]);
%% Reproduce sound
%soundsc(audio_file,Fs_audio);