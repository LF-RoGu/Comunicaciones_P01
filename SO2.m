%% P01_705694
filename= 'Never_Gonna.wav';

fir_order = 100;
Fs_audio = 44100;
Ts_audio = 1/Fs_audio;
t = 10; %t = 10 seconds

samples = [1,t*Fs_audio];

[x,Fs] = audioread(filename,samples);   %obtencion del archivo wav

info = audioinfo(filename);

%% Signal to Mono
xMono = (x(:,1)+ x(:,2))/2;

%% Filtrar la señal 
    %LPF Comunication Channel  
B = 15e3;
Fc = B / (Fs_audio/2);

fc_15KHz = [0 Fc Fc 1];
m = [1 1 0 0];
orden = fir_order;
LPF = fir2(orden, fc_15KHz, m); % LPF (Low Pass Filter)
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
    
    P_noise(i) = (sum(noise.*noise))/numel(noise);
    % Signal to Noise Relation
    SNR_A(i) = P_signal(1)./P_noise(1); 
    % SNR in dB
    SNR_A_dB(i) = 10*log10(SNR_A(i));
end  

    % Add up noise to the signal
xa_noised = xa + noise';
    % Normilized the signal between 1 & -1
xa_noised = xa_noised ./ max(abs(xa_noised));

filename = strcat('AnalogSignal_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
audiowrite(filename, xa_noised, Fs); %obtencion del archivo 

soundsc(xa_noised,Fs_audio);

%% Digital Transmition
