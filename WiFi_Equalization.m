function [ equalized_symbols ] = WiFi_Equalization(rx,H_EST_LT,H_EST_PS,sampUtil,no_OFDM_blocks)
    equalized_symbols = zeros(sampUtil, no_OFDM_blocks);
    for i = 1:no_OFDM_blocks
        H_UTIL = ((no_OFDM_blocks-i)/no_OFDM_blocks).*H_EST_LT + ...
                (i/no_OFDM_blocks).*H_EST_PS;
        equalized_symbols(1:26,i) = rx(1:26,i)./H_UTIL(1:26);
        equalized_symbols(28:end,i) = rx(28:end,i)./H_UTIL(28:end);
    end
end