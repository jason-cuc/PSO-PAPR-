function [PAPR_MPSO,Itermax]= MPSO_PTS( Symbol_ifft2,W,Gn,initial_w,initial_v,c1,c2,Vmax,w,threshold )
%   This Programe compare the min papr based on the PSO technique-based PTS technique
%   according to the article "A Suboptimal PTS Algorithm Based on Particle Swarm 
%   Optimization Technique for PAPR Reduction in OFDM Systems".
%   
%   input:
%       symbol_ifft2: OFDM symbol by partition M subblocks and IFFT
%       W: the set of phase factor
%       Gn: the max iteration number
%       initial_b: initial phase factor matrix
%       initial_v: initial velocity
%       c1,c2: learning factor
%       Vmax: the max velocity
%       w: inertia weight vector
%       delta: decision parameter
%   output:
%       PAPR_PSO: the min PAPR based on PSO_PTS

[M,N] = size(initial_w);                    % the number of particle per generation
w_pbest = zeros(M,N);                       % used to note the particle best position
w_gbest = zeros(M,1);                       % used to note the swarm best position
papr_pbest = inf*ones(1,N);                 % used to note the particle best value
papr_gbest = inf;                           % used to note the swarm best value

current_w = initial_w;                      % update the particle position
current_v = initial_v;                      % update the particle velocity

for ii = 1:1:Gn
    
    Bdata = binary2factor(current_w,W);             % Inverse Mapping
    
    Symbol_ifft = Symbol_ifft2*Bdata;               % compute papr of all particle 
    PowerPerBit = abs(Symbol_ifft).^2;
    PowerMean = mean(PowerPerBit);
    PowerMax  = max(PowerPerBit);
    papr = PowerMax./PowerMean;
    
    update_index1 =  find(papr < papr_pbest) ;      % update the particle best
    w_pbest(:,update_index1) = current_w(:,update_index1);
    papr_pbest(update_index1) = papr(update_index1);
    
    if ( min(papr) < papr_gbest )                   % update the swarm best
        [~,update_index2] = min(papr);   
        w_gbest(:,1) = current_w(:,update_index2);
        papr_gbest = papr(update_index2);
    end
    
    if ( papr_gbest <= threshold )
        break;
    end
    
    next_v = w(ii)*current_v + c1*rand(M,N).*(w_pbest-current_w)...
                                + c2*rand(M,N).*(w_gbest*ones(1,N)-current_w);    % update velocity
    [x_index,y_index] = find( next_v < -Vmax );
    next_v(x_index,y_index) = -Vmax;
    [x_index,y_index] = find( next_v > Vmax );
    next_v(x_index,y_index) = Vmax;
    
    Sv = 1./(1+exp(-next_v));                       % sigma function
    next_w = rand(M,N) < Sv;                        % update position
    
    current_w = next_w;                             % update state
    current_v = next_v;
    
end

PAPR_MPSO = papr_gbest;
Itermax = ii;
    
    


