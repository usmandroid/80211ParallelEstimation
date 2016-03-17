% clear; close all;

% %% PARAMETERS
% Fs = 10e6;                                          % Sampling Frequency
% SNR = 40;                                           % Signal to Noise Ratio
% Pawgn = -50;                                        % Noise Power in dBm
% threshold = 4e-3;                                   % Threshold for signal detection
% channel_model = 'A';                                % Values: A, B, C, D, E
% FO = 20e3;                                          % Frequency Offset (impairment)
% 
% K = 64;                                             % Subcarriers per OFDM Symbol
% sampUtil = 53;                                      % Non-null Subcarriers + DC
% sampXblock = K + 16;                                % Subcarriers + Cyclic Prefic per OFDM Symbol
% nOFDM_symb= 15;                                     % OFDM Symbols per Frame
% 
% %% INPUTS GENERATION
% [ tx_packet, rx_packet , tx_lptot , rx_lptot] = WiFi_inputs();
% 
% tx_preamble1 = tx_lptot(end-63:end,1);
% tx_preamble2 = tx_lptot(end-64-63:end-64,1);
% tx_preamble = (tx_preamble1+tx_preamble2)./2;
% tx_preamble_fft = circshift(fft(tx_preamble,64),26);
% tx_preamble_fft = tx_preamble_fft(1:sampUtil);
% 
% rx_preamble1 = rx_lptot(end-63:end,1);
% rx_preamble2 = rx_lptot(end-64-63:end-64,1);
% rx_preamble = (rx_preamble1+rx_preamble2)./2;
% rx_preamble_fft = circshift(fft(rx_preamble,64),26);
% rx_preamble_fft = rx_preamble_fft(1:sampUtil);
% 
ow2 = sum( (rx_preamble2 - rx_preamble1).*conj(rx_preamble2 - rx_preamble1) )./(2*K);

%% FOURIER MATRIX
F = zeros(K,K);
for f = 1:K
    for t = 1:K
        F(t,f) = exp(-1i*2*pi*(t-1)*(f-1)/K);
    end
end

%% OFDM SYMBOLS EXTRACTION
[ tx_symb, tx_symb_long ] = WiFi_blocks_extraction(tx_packet,K,sampXblock,sampUtil,nOFDM_symb);
[ rx_symb, rx_symb_long ] = WiFi_blocks_extraction(rx_packet,K,sampXblock,sampUtil,nOFDM_symb);

%% CHANNEL ESTIMATION
% LT Linear Square
H_EST_LT_LS = WiFi_channel_estimation_LT_LS(tx_preamble_fft,rx_preamble_fft);

% PS Linear Interpolation
H_EST_PS_Linear = WiFi_channel_estimation_PS_Linear(tx_symb,rx_symb,sampUtil,nOFDM_symb);
% PS Third Order Interpolation
H_EST_PS_Third = WiFi_channel_estimation_PS_Third(tx_symb,rx_symb,sampUtil,nOFDM_symb);
% PS Natural Cubic Interpolation
H_EST_PS_Cubic = WiFi_channel_estimation_PS_Cubic(tx_symb,rx_symb,sampUtil,nOFDM_symb);
% PS Sinc Interpolation
H_EST_PS_Sinc = WiFi_channel_estimation_PS_Sinc(tx_symb,rx_symb,sampUtil,nOFDM_symb);
% PS MMSE
H_EST_PS_MMSE = WiFi_channel_estimation_PS_MMSE(tx_symb,rx_symb,ow2,H_EST_LT_LS,sampUtil,nOFDM_symb);

%% EQUALIZATION
eq_symbols = WiFi_Equalization(rx_symb,H_EST_LT_LS,H_EST_PS_Linear,sampUtil,nOFDM_symb);

%%
tx_symb1 = reshape(tx_symb,[53*15,1]);
rx_symb1 = reshape(rx_symb,[53*15,1]);
