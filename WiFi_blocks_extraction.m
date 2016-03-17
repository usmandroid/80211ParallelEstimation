function [ symb, symb_long ] = WiFi_blocks_extraction(data,K,sampXblock,sampUtil,blocks)
    t_ofdm_blocks = zeros(sampXblock,blocks);
    symb = zeros(sampUtil,blocks);
    symb_long = zeros(K,blocks);
    for nob = 1:blocks
        t_ofdm_blocks(:,nob) = data(1 + sampXblock*(nob-1) : sampXblock*nob);
        symb_long(:,nob) = fft(t_ofdm_blocks(end-63:end,nob));
        symb_long1 = circshift(symb_long(:,nob),26);
        symb(:,nob) = symb_long1(1:sampUtil);
    end
end