%% P01_705694

close all; clc; clear;

%% Transmicion analogica
filename= 'Never_Gonna.wav';
Fs_audio = 44100;
sec = 10;
samples = [1,sec*Fs_audio];
[x,Fs] = audioread(filename,samples);   %obtencion del archivo wav 
info = audioinfo(filename);

%% Pasar a señal mono 
xMono = (x(:,1)+ x(:,2))/2;

%% Filtrar la señal 
%LPF 15 kHz
fc = 15000 / (Fs_audio/2);
B = 15e3;
fc_15KHz = [0 fc fc 1];
m = [1 1 0 0];
orden = 100;
LPF_15KHz = fir2(orden, fc_15KHz, m);
%%fvtool(LPF_15KHz, 'Analysis', 'impulse');
xa = filter(LPF_15KHz, 1, xMono); %señal analogica filtrada
%mean(xa) componente de directa
%min(xa)
%max(xa) %entre -1 y 1 v
%pwelch(xa)

%% Normalizacion de la potencia 
%calcular la potencia
pot = sum(xa.*xa)/numel(xa);
xa = xa/sqrt(pot);
P_signal = var(xa); %verificar la potencia 

%% SNR
B =15e3 ; % Ancho de Banda del filtro receptor
N0 = 1./(B.*10.^(0:0.3:3)); % Vector PSD del ruido
P_noise = B*N0; % Vector de Pot del Ruido Filtrado
P_noise_dB = 10.*log10(P_noise); % Pot. Ruido en Decibeles
SNR_A = P_signal ./ P_noise; % Relacion Señal a Ruido
SNR_A_dB = 10*log10(SNR_A); % SNR en dB

%ruido AWGN 
for i=1:length(N0)
noise = sqrt(P_noise(i)) * randn(size(xa)); %Muestras de ruido
xa_noised = xa + noise; 
xa_noised = xa_noised ./ max(abs(xa_noised));%normalizar entre 1 y -1

filename = strcat('AnalogSignal_P1_', num2str(i), '_', num2str(round(SNR_A_dB(i)), 4), 'dB','.wav');
audiowrite(filename, xa_noised, Fs); %obtencion del archivo 
end 
soundsc(xa_noised,Fs_audio);

%% señal digital 
 %Source: https://www.mathworks.com/matlabcentral/answers/75379-how-to-convert-an-input-sine-wave-into-an-8-bit-digital-signal
    %Source: https://www.mathworks.com/matlabcentral/answers/103315-how-to-convert-digital-data-into-analog-data-using-matlab-code
    % Convert values from xa_noised to bin

M = 2^info_bps - 1;

level = max(xa_noised);
level = level/M;

step = level;

data = abs(xa_noised);

xa2bi = fix(data/step);
    %For Default it is the 'right-msb'
b = de2bi(xa2bi,info_bps,'left-msb');
    %Transponse operation
b = b';

