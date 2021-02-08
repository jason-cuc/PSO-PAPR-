function IQvalue=get80216map(Qam) 
% This function is used to calculate constellation IQ vector 
% Qam=4 -- QPSK, Qam=16 -- 16-QAM, Qam=64 -- 64QAM 
switch Qam, 
    case 4, 
        IQvalue = [1+1i 1-1i -1+1i -1-1i]/sqrt(2); 
    case 16, 
        col1 = [1+1i 1+3i 1-1i 1-3i];  % column 1 of constellation 
        IQvalue = [col1 col1+2]; 
        IQvalue = [IQvalue conj(-IQvalue)]/sqrt(10); 
    case 64, 
        quad1 = [3+3i 3+1i 3+5i 3+7i]; % quadrant 1 
        col1 = [quad1 conj(quad1)]; 
        IQvalue = [col1 col1-2]; 
        IQvalue = [IQvalue col1+2]; 
        IQvalue = [IQvalue col1+4]; 
        IQvalue = [IQvalue conj(-IQvalue)]/sqrt(42); 
end