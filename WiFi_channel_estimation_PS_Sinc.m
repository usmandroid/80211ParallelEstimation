function [ H_EST ] = WiFi_channel_estimation_PS_Sinc(tx,rx,sampUtil,no_OFDM_blocks)
    H_EST_ML_PILOT_SINC = zeros(sampUtil,4);
    sinc1 = zeros(sampUtil,1);
    sinc2 = zeros(sampUtil,1);
    sinc3 = zeros(sampUtil,1);
    sinc4 = zeros(sampUtil,1);
    for i = 1:no_OFDM_blocks
        rx_pilots = [rx(6,i);rx(20,i);rx(34,i);rx(48,i)];
        tx_pilots = [tx(6,i);tx(20,i);tx(34,i);tx(48,i)];
        H_PILOTS = rx_pilots./tx_pilots;
        for k = 1:sampUtil
            sinc1(k) = H_PILOTS(1)*sinc((k-6)/(20-6));
            sinc2(k) = H_PILOTS(2)*sinc((k-20)/(20-6));
            sinc3(k) = H_PILOTS(3)*sinc((k-34)/(34-20));
            sinc4(k) = H_PILOTS(4)*sinc((k-48)/(48-34));
        end
        H_EST_ML_PILOT_SINC(:,i) = sinc1+sinc2+sinc3+sinc4;
    end
    H_EST = ( H_EST_ML_PILOT_SINC(:,1)+H_EST_ML_PILOT_SINC(:,2)+H_EST_ML_PILOT_SINC(:,3)+H_EST_ML_PILOT_SINC(:,4) )./4;
end