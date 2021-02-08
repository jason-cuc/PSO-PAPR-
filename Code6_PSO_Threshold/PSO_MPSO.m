%   This Programe compares the performances using the PSO technique-based PTS techniques
%   while different generation number.
%   
%   Author: dtkong
%   Version: 1.0
%   Time: 2011.03.29

close all;
clear all;
clc;

%============================= setting common parameter =================================%
NumCarr = 1024;                 % the number of transmission subcarriers
NumSymb = 1e4;                  % the number of symbols
mapsize = 2;                    % using QPSK modulation
V = 16;                         % the number of subblocks 
OverSampleRate = 4;             % over sample rate
Partition = 1;               	% the way of partition:1 -> adjacency partition;2 -> interlaced partition;3 -> random partition
W = 1;                          % the log of the length of weighting factor set
%========================================================================================%

%============================ setting intial parameter for PSO_PTS ======================%
Num_Particle = 10;              % the number of particles per generation
Gn = 10;                        % the max iteration number
threshold = 6.7;                % the threshold value of PAPR
c1 = 2; c2 = 2;                 % learning factor
Vmax = 0.2;                     % the max velocity
wmax = 0.9;                     % the max inertia weight vector
wmin = 0.4;                     % the min inertia weight vector
w = wmax-(wmax-wmin)/Gn*(1:Gn); % inertia weight vector
initial_w = randi([0,1],[W*V,Num_Particle]);                % initial position  
v_min = -Vmax;  v_max = Vmax;
initial_v = v_min + (v_max-v_min).*rand(W*V,Num_Particle);  % inital volcity
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
PAPROnIPTS = zeros(1,NumSymb);              % use to note the papr with IPTS
PAPR_PSO = zeros(1,NumSymb);                % use to note the papr with PSO_PTS while PSO
PAPR_MPSO = zeros(1,NumSymb);               % use to note the papr with PSO_PTS while MPSO
PAPR_Itermax = zeros(1,NumSymb);            % use to the note the number of Iteration while with threshold

PAPR0 = 0:0.25:12;                          % the vaule range of PAPR0

CntNoPTS = zeros(1,length(PAPR0));          % use to note the value of CCDF without PTS
Cnt_IPTS = zeros(1,length(PAPR0));          % use to note the value of CCDF with IPTS
Cnt_PSO = zeros(1,length(PAPR0));           % use to note the value of CCDF with PSO
Cnt_MPSO = zeros(1,length(PAPR0));          % use to note the value of CCDF with MPSO

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

  	%------------------------------- Partition OFDM Symbol ------------------------------%
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
  	%------------------------------------------------------------------------------------%
    
    %---------------------------------- IPTS --------------------------------------------%
    PAPROnIPTS(1,n) = IPTS( Symbol_ifft2,V,2^W,[-1 1] );
    %------------------------------------------------------------------------------------%


  	%--------------------- PAPR with PSO_PTS while PSO and MPSO -------------------------%
    PAPR_PSO(1,n) = PSO_PTS( Symbol_ifft2,W,Gn,initial_w,initial_v,c1,c2,Vmax,w );
    [PAPR_MPSO(1,n),PAPR_Itermax(1,n)] = MPSO_PTS( Symbol_ifft2,W,Gn,initial_w,initial_v,c1,c2,Vmax,w,threshold );
    %--------------------------------------------------------------------------------%
    n
end
%----------------------------------------------------------------------------------------%

PAPR_Iter = sum(PAPR_Itermax)/NumSymb

PAPRNoPTS = 10*log10(PAPRNoPTS);
PAPROnIPTS = 10*log10(PAPROnIPTS);
PAPR_PSO = 10*log10(PAPR_PSO);
PAPR_MPSO = 10*log10(PAPR_MPSO);

for k = 1:1:length(PAPR0)
    CntNoPTS(k) = sum( PAPRNoPTS > PAPR0(k) );
    Cnt_IPTS(k) = sum( PAPROnIPTS > PAPR0(k) );
    Cnt_PSO(k) = sum( PAPR_PSO > PAPR0(k) );
    Cnt_MPSO(k) = sum( PAPR_MPSO > PAPR0(k) );
end

CntNoPTS = CntNoPTS./NumSymb;
Cnt_IPTS = Cnt_IPTS./NumSymb;
Cnt_PSO = Cnt_PSO./NumSymb;
Cnt_MPSO = Cnt_MPSO./NumSymb;

figure(1)
semilogy(PAPR0,CntNoPTS,'-k.');
hold on
semilogy(PAPR0,Cnt_IPTS,'-go');
hold on
semilogy(PAPR0,Cnt_PSO,'-r*');
hold on
semilogy(PAPR0,Cnt_MPSO,'-b+');
grid on;
legend('without PTS','IPTS','PSO','PSO with threshold');
xlabel('PAPR0(dB)'),ylabel('Pr(PAPR>PAPR0)');
title(['Comparision of PSO-PTS with PSO and PSO with threshold']);
xlim([6,14]);
hold off



