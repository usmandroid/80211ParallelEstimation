function [ H_EST ] = WiFi_channel_estimation_PS_Cubic(tx,rx,sampUtil,no_OFDM_blocks)
    H_EST_ML_PILOT_THIRD = zeros(sampUtil,no_OFDM_blocks);
    for i = 1:no_OFDM_blocks
        rx_pilots = [rx(6,i);rx(20,i);rx(34,i);rx(48,i)];
        tx_pilots = [tx(6,i);tx(20,i);tx(34,i);tx(48,i)];
        H_PILOTS = rx_pilots./tx_pilots;
        f0    = H_PILOTS(1);
        f01   = (H_PILOTS(2)-H_PILOTS(1))/(20-6);
        f12   = (H_PILOTS(3)-H_PILOTS(2))/(34-20);
        f23   = (H_PILOTS(4)-H_PILOTS(3))/(48-34);
        f012  = (f12-f01)/(34-6);
        f123  = (f23-f12)/(48-20);
        f0123 = (f123-f012)/(48-6);
        for k = 1:sampUtil
            H_EST_ML_PILOT_THIRD(k,i) = f0 + f01*(k-6) + f012*(k-6)*(k-20) + f0123*(k-6)*(k-20)*(k-34);
        end
    end
    H_EST = ( H_EST_ML_PILOT_THIRD(:,1)+H_EST_ML_PILOT_THIRD(:,2)+H_EST_ML_PILOT_THIRD(:,3)+H_EST_ML_PILOT_THIRD(:,4) )./4;
end