function [ H_EST ] = WiFi_channel_estimation_PS_Linear(tx,rx,sampUtil,no_OFDM_blocks)
    H_EST = zeros(sampUtil,4);
    for i = 1:no_OFDM_blocks
        rx_pilots = [rx(6,i);rx(20,i);rx(34,i);rx(48,i)];
        tx_pilots = [tx(6,i);tx(20,i);tx(34,i);tx(48,i)];
        H_PILOTS = rx_pilots./tx_pilots;
        for k = 1:sampUtil
            if(k<6)
                alpha = (k-6)/(20-6);
                H_EST(k,i) = H_PILOTS(1)+( (H_PILOTS(2)-H_PILOTS(1) )*alpha );
            elseif ((k>=6) && (k<20))
                alpha = (k-6)/(20-6);
                H_EST(k,i) = H_PILOTS(1)+( (H_PILOTS(2)-H_PILOTS(1) )*alpha );
            elseif ((k>=20) && (k<34))
                alpha = (k-20)/(34-20);
                H_EST(k,i) = H_PILOTS(2)+( (H_PILOTS(3)-H_PILOTS(2) )*alpha );
            elseif ((k>=34) && (k<48))
                alpha = (k-34)/(48-34);
                H_EST(k,i) = H_PILOTS(3)+( (H_PILOTS(4)-H_PILOTS(3) )*alpha );
            elseif (k>=48)
                alpha = (k-34)/(48-34);
                H_EST(k,i) = H_PILOTS(3)+( (H_PILOTS(4)-H_PILOTS(3) )*alpha );
            end
        end
    end
    H_EST = ( H_EST(:,1)+H_EST(:,2)+H_EST(:,3)+H_EST(:,4) )./4;
end