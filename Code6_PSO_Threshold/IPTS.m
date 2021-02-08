function PAPR_IPTS = IPTS( Symbol_ifft2,M,W,weight_factor )
%   This Programe compares the papr of ofdm system using the iterative algorithm according
%   the article "Peak-to-average Power Ratio reduction of an OFDM signal using partial 
%   transmit sequences".
%
%   input:
%       symbol_ifft2: OFDM symbol by partition M subblocks and IFFT
%       M: the number of subblocks
%       W: the element number of phase factor set
%       weight_factor: the element set of phase rotation factors
%   output:
%       PAPR_IPTS: the minimun PAPR according this method

bdata = ones(M,1);              % the phase rotation factors vector

for ii = 1:1:M
    Bdata = bdata*ones(1,W);            % the all combinations
    Bdata(ii,:) = weight_factor;
    Symbol_ifft = Symbol_ifft2*Bdata;   % output OFDM symbol
    PowerPerBit = abs(Symbol_ifft).^2;
    PowerMean = mean(PowerPerBit);
    PowerMax  = max(PowerPerBit);
    PAPR = PowerMax./PowerMean;         % compute papr for different combinations
    [papr_min,min_index] = min(PAPR);
    bdata = Bdata(:,min_index);
end

PAPR_IPTS = papr_min;