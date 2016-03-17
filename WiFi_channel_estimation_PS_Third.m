function [ H_EST ] = WiFi_channel_estimation_PS_Third(tx,rx,sampUtil,no_OFDM_blocks)
    H_EST_ML_PILOT_SPLINE = zeros(sampUtil,no_OFDM_blocks);
    for i = 1:no_OFDM_blocks
        rx_pilots = [rx(6,i);rx(20,i);rx(34,i);rx(48,i)];
        tx_pilots = [tx(6,i);tx(20,i);tx(34,i);tx(48,i)];
        H_PILOTS = rx_pilots./tx_pilots;
        xx = linspace(1,sampUtil,sampUtil);
        x = [6 20 34 48];
        y = [H_PILOTS(1) H_PILOTS(2) H_PILOTS(3) H_PILOTS(4)];
        H_EST_ML_PILOT_SPLINE(:,i) = csapi(x,y,xx);
    end
    H_EST = ( H_EST_ML_PILOT_SPLINE(:,1)+H_EST_ML_PILOT_SPLINE(:,2)+H_EST_ML_PILOT_SPLINE(:,3)+H_EST_ML_PILOT_SPLINE(:,4) )./4;
end