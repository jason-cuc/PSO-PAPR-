%   This Programe compares the performances using the PSO technique-based PTS techniques
%   while different generation number.
%   
%   Author: dtkong
%   Version: 1.0
%   Time: 2010.12.1

close all;
clear all;
clc;

%============================= setting common parameter =================================%
NumCarr = 1024;                 % the number of transmission subcarriers
NumSymb = 1e5;                  % the number of symbols
mapsize = 2;                    % using QPSK modulation
VV = [2 4 8 16];                % the number of subblocks 
OverSampleRate = 4;             % over sample rate
Partition = 1;                  % the way of parition:1 -> adjacency partition;2 -> interlaced partition
W = 1;                          % the log of the length of weighting factor set
%========================================================================================%

%============================ setting intial parameter for PSO_PTS ======================%
Num_Particle = 10;              % the number of particles per generation
Gn = 30;                        % the max iteration number
c1 = 2; c2 = 2;                 % learning factor
Vmax = 0.2;                     % the max velocity
wmax = 0.9;                     % the max inertia weight vector
wmin = 0.4;                     % the min inertia weight vector
w = wmax-(wmax-wmin)/Gn*(1:Gn); % inertia weight vector
v_min = -Vmax;  v_max = Vmax;
%========================================================================================%

%========================== initialization data of system ===============================%
MapSymb = get80216map( 2^mapsize );         % initial mapping method
Initial_Pattern = ones( NumCarr,1 );        % initial all of position as 1
Initial_Pattern([1:28,end-26:end]) = 0;     % initial guard band as 0
Initial_Pattern(45:24:end-27) = 4;          % initial pilot position as 4
Position_pilot = find(Initial_Pattern==4);  % the position of pilot
Position_data  = find(Initial_Pattern==1);  % the position of transimitted data
Length_pilot = length(Position_pilot);      % the length of pilot
Length_data  = length(Position_data);       % the length of transimitted data

PAPRNoPTS = zeros(1,NumSymb);               % use to note the papr without PTS
PAPR_PSO = zeros(length(VV),NumSymb);    	% use to note the papr with PSO_PTS

PAPR0 = 0:0.25:12;                          % the vaule range of PAPR0

CntNoPTS = zeros(1,length(PAPR0));          % use to note the value of CCDF without PTS
Cnt_PSO = zeros(length(VV),length(PAPR0)); 	% use to note the value of CCDF with PSO_PTS
%========================================================================================%

%---------------------------------- Computing PAPR  -------------------------------------%
for n = 1:1:NumSymb
    
    Datatx = floor(rand(Length_data,1)*(2^mapsize));
    DataMap = MapSymb(Datatx+1);
    Symbol_tx = Initial_Pattern;
    Symbol_tx(Position_pilot) = round(rand(Length_pilot,1));
    Symbol_tx(Position_data)  = DataMap;
    
    %---------------------------------- PAPR without PTS --------------------------------%
    Symbol_ifft = ifft([Symbol_tx(1:NumCarr/2);zeros(NumCarr*(OverSampleRate-1),1);...
                                                            Symbol_tx(NumCarr/2+1:end)]);
    PowerPerBit = abs(Symbol_ifft).^2;
    PowerMean = mean(PowerPerBit);
    PowerMax  = max(PowerPerBit);
    
    PAPRNoPTS(1,n) = PowerMax/PowerMean;
    %------------------------------------------------------------------------------------%
    
    for kk = 1:1:length(VV)
        %------------------------------- Partition OFDM Symbol --------------------------%
        V = VV(kk);
        Symbol_block = zeros( NumCarr,V );
        Symbol_ifft2 = zeros( NumCarr*OverSampleRate,V );
        for v = 1:1:V
            if ( Partition == 1 )
                Symbol_block(((v-1)*NumCarr/V+1):(v*NumCarr/V),v) = Symbol_tx(((v-1)*NumCarr/V+1):(v*NumCarr/V));
            elseif ( Partition == 2 )
                Symbol_block(v:V:NumCarr,v) = Symbol_tx(v:V:NumCarr);
            end
            Symbol_ifft2(:,v) = ifft([Symbol_block(1:NumCarr/2,v);zeros(NumCarr*(OverSampleRate-1),1);...
                                                                            Symbol_block(NumCarr/2+1:end,v)]);
        end
        %--------------------------------------------------------------------------------%
    
        %----------------------------- PAPR with PSO_PTS --------------------------------%
        initial_w = randi([0,1],[W*V,Num_Particle]);                % initial position  
        initial_v = v_min + (v_max-v_min).*rand(W*V,Num_Particle);  % inital volcity
        PAPR_PSO(kk,n) = PSO_PTS( Symbol_ifft2,W,Gn,initial_w,initial_v,c1,c2,Vmax,w );
        %--------------------------------------------------------------------------------%
    end

    n
end
%----------------------------------------------------------------------------------------%

PAPRNoPTS = 10*log10(PAPRNoPTS);
PAPR_PSO = 10*log10(PAPR_PSO);

for k = 1:1:length(PAPR0)
    CntNoPTS(k) = sum( PAPRNoPTS > PAPR0(k) );
    for ii = 1:1:length(VV)
        Cnt_PSO(ii,k) = sum( PAPR_PSO(ii,:) > PAPR0(k) );
    end
end

CntNoPTS = CntNoPTS./NumSymb;
Cnt_PSO = Cnt_PSO./NumSymb;

figure(1)
semilogy(PAPR0,CntNoPTS,'-r*');
hold on
semilogy(PAPR0,Cnt_PSO(1,:),'-ks');
hold on;
semilogy(PAPR0,Cnt_PSO(2,:),'-cd');
hold on;
semilogy(PAPR0,Cnt_PSO(3,:),'-rh');
hold on;
semilogy(PAPR0,Cnt_PSO(4,:),'-mp');
grid on;
legend('without PTS','V = 2','V = 4','V = 8','V = 16');
xlabel('PAPR0(dB)'),ylabel('Pr(PAPR>PAPR0)');
title('Comparision of based PSO PTS with different V')
xlim([5 12]);
hold off



