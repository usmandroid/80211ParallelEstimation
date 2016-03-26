function [ H_EST_MMSE_PILOT_AV ] = WiFi_channel_estimation_PS_MMSE(tx,rx,ow2,H_EST,sampUtil,no_OFDM_blocks)
    H_EST_MMSE_PILOT = zeros(sampUtil,4);
    tx(1:5,:)=0;
    tx(7:19,:)=0;
    tx(21:33,:)=0;
    tx(35:47,:)=0;
    tx(47:end,:)=0;
    rx(1:5,:)=0;
    rx(7:19,:)=0;
    rx(21:33,:)=0;
    rx(35:47,:)=0;
    rx(47:end,:)=0;

    F = zeros(sampUtil,sampUtil);
    for f = 1:sampUtil
        for t = 1:sampUtil
            F(t,f) = exp(-1i*2*pi*(t-1)*(f-1)/sampUtil);
        end
    end

    Rhh = ifft(H_EST,sampUtil)*ifft(H_EST,sampUtil)';

    for i = 1:no_OFDM_blocks
        X4 = diag(tx(:,i));
        Rhy = Rhh*F'*X4;
        Ryy = X4*F*Rhh*F'*X4' + ow2*eye(sampUtil);
        H_EST_MMSE_PILOT(:,i) = F*Rhy*pinv(Ryy)*rx(:,i);
    end
    H_EST_MMSE_PILOT_AV = ( H_EST_MMSE_PILOT(:,1)+H_EST_MMSE_PILOT(:,2)+H_EST_MMSE_PILOT(:,3)+H_EST_MMSE_PILOT(:,4) )./4;
end