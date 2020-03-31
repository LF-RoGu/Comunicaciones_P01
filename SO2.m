%% P01_705694_710786

close all; clc; clear;

filename= 'Never_Gonna.wav';    %audio to be reproduce

fir_order = 100;    %filter order (15kHz) 
Fs_audio = 44100;   %sample frequency of audio
Ts_audio = 1/Fs_audio;  %period of the sample audio
t = 10; %t = 10 seconds
samples = [1,t*Fs_audio]; %vector of samples
[x,Fs] = audioread(filename,samples);   %get wav file
info = audioinfo(filename); %get info of the wav file
info_bps = info.BitsPerSample; %info_bps (info_bits per sample)
    
%Signal to Mono
xMono = (x(:,1)+ x(:,2))/2; %convert audio stero to audio mono

%% Filtrar la señal 
    %LPF Comunication Channel  
B = 15e3; %band pass of filter
Fc = B / (Fs_audio/2); %Normalized Freq 0.6803

fc_15KHz = [0 Fc Fc 1]; %Freq of the filter
m = [1 1 0 0];  %Magnitud of filter
order = fir_order; %order of teh filter
LPF = fir2(order, fc_15KHz, m); % LPF (Low Pass Filter)
xa = filter(LPF, 1, xMono); %analog signal filtered

%% Normalize the Power
pot = sum(xa.*xa)/numel(xa);  % Calculate the power of analog signal
xa = xa/sqrt(pot);  %Unitari normalization of power
P_signal = var(xa); %Check power
%% SNR
N0 = 1./(B.*10.^(0:0.3:3)); % PSD vector of noise
P_noise = B*N0; % Power of filtered noise
P_noise_dB = 10.*log10(P_noise); % Power of noise in dB 
SNR_A = P_signal ./ P_noise; % Signal to Noise ratio 
SNR_A_dB = 10*log10(SNR_A); % SNR in dB

%% AWGN (White noise all frequency)
% Add up noise to the analog signal, for each value of power density of noise
for i = 1 : numel(N0)
    noise = sqrt(P_noise(i)) * randn(1, numel(xa)); % noise samples for each value noise power
    xa_noised = xa + noise';     % Add up noise to the signal  
    xa_noised = xa_noised ./ max(abs(xa_noised)); % Normilized the signal between 1 & -1
    %save name of the new file that contain the audio signal + noise
    filename = strcat('AnalogSignal_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
    audiowrite(filename, xa_noised, Fs); %get new file
end 

% Reproduce sound
%soundsc(xMono,Fs_audio);

%% Digital Transmition
    %Source: https://www.mathworks.com/matlabcentral/answers/75379-how-to-convert-an-input-sine-wave-into-an-8-bit-digital-signal
    %Source: https://www.mathworks.com/matlabcentral/answers/103315-how-to-convert-digital-data-into-analog-data-using-matlab-code
    % Convert values from xMono to bin

swing = (2^info_bps-1)/2; %offset to be added
xMono2bi = round(xMono*swing+swing);    %add offset to digital audio
    %For Default it is the 'right-msb'
b = de2bi(xMono2bi,info_bps,'left-msb');    %pass decimal to binari    
b = b'; %Transponse operation
bits_tx = b(:); %Concatenate the rest of the bits, to make it a single vector

%This goes if beta it is equal to ZERO
%B = 15e3;
%Rb_max = B2*2;
%% Shapping Pulse
beta = 0.35; 
D = 6; 
Fs = 96e3; 
Rb = (2*B) / (1+beta); 
Ts = 1/Fs; 
mp = ceil(Fs/Rb); %mp = 5;
Tp = 1/(Fs/mp); 
energy = Tp; 

[p,t] = rcpulse(beta, D, Tp, Ts, 'srrc', energy); %Genetare
%stem(p);
p_energy = Ts*sum(p.*p);    %Calculate the power of the SRRC pulse
p = p/sqrt(p_energy); %Unitari Normalization of SRRC pulse power
%% Code Line
    %Polar Signal NRZ
PBP = p; % PSPB (Polar Base Pulse)
PS_s = bits_tx; %auxiliar vector that contain bits information
PS_s( PS_s == 0 ) = -1; %convert binary 0 to -1 v
s = zeros(1,numel(PS_s)*mp);    %Vector with the total size
s(1:mp:end) = PS_s;     %Impulse train
PSLC = conv(PBP, s) / mp; %PSLC (Polar Signal Line Code), Pulse train 

%plot(PSLC(1:mp*20));

    %Normilize pulse
PSLC = PSLC/sqrt( sum(PSLC.*PSLC)/numel(PSLC) ); %Unitari Normalization of Line code power
P_PSLC_d = var(PSLC); %Power of Line code digital




%% SNR Digital 
N0 = 1./(B.*10.^(0:0.3:3)); % PSD vector of noise
P_noise_d = B*N0; % Power of filtered noise
P_noise_dB_d = 10.*log10(P_noise_d); % Power of noise in dB 
SNR_D = P_PSLC_d ./ P_noise_d; % Signal to Noise ratio
SNR_D_dB = 10*log10(SNR_D); % SNR in dB

%% Match Filter
%Match filter
%For the match filter we used the same form of the base pulse (SRRC)
[h_srrc t] = rcpulse(beta,D,Tp,Ts,'srrc',energy); %Response of the match Filter
hsrrc_energy = Ts*sum(h_srrc.*h_srrc); %Calculate the energy of the filter
h_srrc = h_srrc/sqrt(hsrrc_energy); %Normalized power of the filter

%% Add noise to Line Code 
% AWGN (White noise all frequency) 
%Generate noise vector and add to polar NRZ line code with SRRC as pulse
%base

%We need to add up noise, filtered, and recover the information inside a for loop because each
%iteration represent one noise power, so is diferent eveery time. 
for i = 1 : numel(N0) 
    
        %noise samples
    noise = sqrt(P_noise(i)) * randn(1, numel(PSLC)); %created noise vector
    
    PSLC_noised = PSLC + noise; % Add up noise to the signal
    
    % filtered line code + noise with the 15kh LPF

        %Conv function between the filter and the signal with AWGN
    PRS_rx = conv(LPF,PSLC_noised); %PSRS (Polar Received Signal)
   %figure; pwelch(PRS_rx,500,300,500,Fs,'power'); title('Polar Fc 0.6803');

    filename_t = strcat('PSLC', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
    
    figure; 
    stem(PSLC_noised(1:mp*10));
    hold on;
    plot(PSLC(1:mp*20));title(filename_t); 
    hold on; 
    plot(PRS_rx(1:mp*20));
    legend('PSLC_noised','PSLC','PRS_rx');

        %filtered the recieved line code in the match filter
    match_filter_recived = conv(h_srrc,PRS_rx);
    
    %Wake the sampling and desition of the line code recived.
    %We have order/2 for the filter delay.
    %We have mp/2 for the mp size, so we can get the signal in the middle
    %of the signal.
    %Take a sample each mp/2 in time domain to recover the data.
    rx_match_filter = match_filter_recived( ( mp/2 + ((order/2)+(numel(p)/2)) : mp : (end - ((order/2)+(numel(p)/2)))));

%en este punto, se obtiene la información recivada y pasada por el filtro
%que emula el canal y el match filter (recepcion)  
    audio_rx = sign(rx_match_filter); %Obtenemos valores entre 1 y -1 del arreglo

    bits_rx = (audio_rx+1)/2; %Normalizamos para que los valores sean 1 o 0
        %Por error de acoplamiento de la señal original, eliminamos las
        %primeras 4 muestras y las 3 ultimas.
        %Antes de esto teniamos error de 50% (por estadistica, pues esta
        %bien y mal), despues de arreglar esto, el error fue del 8%
    bits_rx = bits_rx(4:end - 3);
    bits_rx = bits_rx';
        %Checksum de la señal original y la recuperada, obtenemos error del
        %8%
    error = (sum( xor(bits_tx,bits_rx) ) / numel(bits_rx))*100;

%Convert vector to mat
    mat = vec2mat(bits_rx,16);
%Convert mat to dec
    audio_bin2dec = bi2de(mat,'left-msb');

    audio_mat = vec2mat(audio_bin2dec, 1 , Fs_audio);

%Normilize audio between 1 & -1
    %audio_mat = 2*mat2gray(audio_mat) - 1;
    audio_file = normalize(audio_mat, 'range', [-1 1]);

    filename = strcat('DigitalSignalRx_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
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

EP = comm.EyeDiagram('SampleRate',Fs*mp,'SamplesPerSymbol',mp);

PRS_rx = PRS_rx';
%Eye Pattern Received Signal
EP(PRS_rx);


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
%% Ejercicio 2 
% SNR Digital 
B_2 = 476280; 
%N0_2 = 1./(B_2.*10.^(0:0.3:3)); % Vector PSD del ruido
P_noise_d_2 = B_2*N0; % Vector de Pot del Ruido Filtrado
P_noise_dB_d_2 = 10.*log10(P_noise_d_2); % Pot. Ruido en Decibeles
SNR_D_2 = P_PSLC_d ./ P_noise_d_2; % Relacion Señal a Ruido
SNR_D_dB_2 = 10*log10(SNR_D_2); % SNR en dB

%% Ejercicio 3 
% comprobacion del n0 
SNR_D_3 = 1/0.03162;
SNR_D3_dB_3 = 10*log10(SNR_D_3);