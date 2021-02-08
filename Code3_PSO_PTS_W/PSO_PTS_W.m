%   This Programe compares the performances using the PSO technique-based PTS techniques
%   while different phase factor set number.
%   
%   Author: dtkong
%   Version: 1.0
%   Time: 2011.03.29

close all;
clear all;
clc;

%============================= setting common parameter =================================%
NumCarr = 1024;                 % the number of transmission subcarriers
NumSymb = 1e5;                  % the number of symbols
mapsize = 2;                    % using QPSK modulation
V = 4;                          % the number of subblocks 
OverSampleRate = 4;             % over sample rate
Partition = 1;                  % the way of partition:1 -> adjacency partition;2 -> interlaced partition;3 -> random partition
WW = [1,2,3];                   % the log of the length of weighting factor set
%========================================================================================%

%============================ setting intial parameter for PSO_PTS ======================%
Num_Particle = 10;              % the number of particles per generation
Gn = 20;                        % the max iteration number
c1 = 2; c2 = 2;                 % learning factor
Vmax = 0.2;                    	% the max velocity
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

% PAPRNoPTS = zeros(1,NumSymb);               % use to note the papr without PTS
% PAPROnPSO_W = zeros(length(WW),NumSymb);    % use to note the papr with PSO_PTS while different W

load PAPRNoPTS.mat;
load PAPROnPSO_W.mat;

PAPR0 = 0:0.25:12;                              % the vaule range of PAPR0

CntNoPTS = zeros(1,length(PAPR0));              % use to note the value of CCDF without PTS
CntOnPSO_W = zeros(length(WW),length(PAPR0));   % use to note the value of CCDF with different W

%========================================================================================%

%---------------------------------- Computing PAPR  -------------------------------------%
for kk = 3:1:length(WW)
    
    W = WW(kk);                                                 % the number of transmission subcarriers
    initial_w = randi([0,1],[W*V,Num_Particle]);                % initial position  
    initial_v = v_min + (v_max-v_min).*rand(W*V,Num_Particle);  % inital volcity
    
    if ( kk == 3 )
        n0 = 10000;
    else
        error('error');
    end
    
    for n = n0:1:NumSymb

        Datatx = floor(rand(Length_data,1)*(2^mapsize));
        DataMap = MapSymb(Datatx+1);
        Symbol_tx = Initial_Pattern;
        Symbol_tx(Position_pilot) = round(rand(Length_pilot,1));
        Symbol_tx(Position_data)  = DataMap;

        if ( kk == 1 )
            %------------------------------ PAPR without PTS ----------------------------%
            Symbol_ifft = ifft([Symbol_tx(1:NumCarr/2);zeros(NumCarr*(OverSampleRate-1),1);...
                                                                    Symbol_tx(NumCarr/2+1:end)]);
            PowerPerBit = abs(Symbol_ifft).^2;
            PowerMean = mean(PowerPerBit);
            PowerMax  = max(PowerPerBit);

            PAPRNoPTS(1,n) = PowerMax/PowerMean;
            %----------------------------------------------------------------------------%
        end

        %------------------------------- Partition OFDM Symbol --------------------------%
        Symbol_block = zeros( NumCarr,V );
        Symbol_ifft2 = zeros( NumCarr*OverSampleRate,V );
        Index = randperm( NumCarr );
        for v = 1:1:V
            if ( Partition == 1 )
                Symbol_block(((v-1)*NumCarr/V+1):(v*NumCarr/V),v) = Symbol_tx(((v-1)*NumCarr/V+1):(v*NumCarr/V));
            elseif ( Partition == 2 )
                Symbol_block(v:V:NumCarr,v) = Symbol_tx(v:V:NumCarr);
            elseif ( Partition == 3 )
                Symbol_block(Index(v:V:NumCarr),v) = Symbol_tx(Index(v:V:NumCarr));
            end
            Symbol_ifft2(:,v) = ifft([Symbol_block(1:NumCarr/2,v);zeros(NumCarr*(OverSampleRate-1),1);...
                                                            Symbol_block(NumCarr/2+1:end,v)]);
        end
        %--------------------------------------------------------------------------------%

        %---------------- PAPR with PSO_PTS while different Paritition ------------------%
        PAPROnPSO_W(kk,n) = PSO_PTS( Symbol_ifft2,W,Gn,initial_w,initial_v,c1,c2,Vmax,w ); 
        %--------------------------------------------------------------------------------%
        [kk,n]

    end
end
%----------------------------------------------------------------------------------------%

PAPRNoPTS = 10*log10(PAPRNoPTS);
PAPROnPSO_W = 10*log10(PAPROnPSO_W);

for k = 1:1:length(PAPR0)
    CntNoPTS(k) = sum( PAPRNoPTS > PAPR0(k) );
    for ii = 1:1:length(WW)
        CntOnPSO_W(ii,k) = sum( PAPROnPSO_W(ii,:) > PAPR0(k) );
    end
end

CntNoPTS = CntNoPTS./NumSymb;
CntOnPSO_W = CntOnPSO_W./NumSymb;

figure(1)
semilogy(PAPR0,CntNoPTS,'-k.');
hold on
semilogy(PAPR0,CntOnPSO_W(1,:),'-r*');
hold on
semilogy(PAPR0,CntOnPSO_W(2,:),'-b+');
hold on
semilogy(PAPR0,CntOnPSO_W(3,:),'-kx');
grid on;
legend('without PTS','W = 2','W = 4','W = 8');
xlabel('PAPR0(dB)'),ylabel('Pr(PAPR>PAPR0)');
title(['Comparision of PSO-PTS while different phase factor set']);
xlim([4,14]);
hold off



