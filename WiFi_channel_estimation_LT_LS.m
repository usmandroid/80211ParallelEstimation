function [ H_EST ] = WiFi_channel_estimation_LT_LS(fft_tx,fft_rx)
    H_EST_ML1 = ( conj(fft_tx(1:26)).*fft_rx(1:26) ) ./ ( conj(fft_tx(1:26)).*fft_tx(1:26) );
    H_EST_ML2 = ( conj(fft_tx(28:end)).*fft_rx(28:end) ) ./ ( conj(fft_tx(28:end)).*fft_tx(28:end) );
    H_EST = [H_EST_ML1;0;H_EST_ML2];
end