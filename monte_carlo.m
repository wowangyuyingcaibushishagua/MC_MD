
function [cdMoved] = monte_carlo(phys_cond,phys_const,index,cdInit,charge,distance_step,L)
% Cdinit stands for initial coordinates, numParticle stands for
% numbers of particle.
cdMoved = cdInit;
% numerator of the formula calculating potential
charge_index = charge(index);
charge(index) = [];
chargemulti = charge_index*charge;
% generate reasonable move
while true
    % generate random vector of movement
    delta = distance_step * normrnd(0, 1 , [1 , 3]);
    cdMoved(index,:) = cdMoved(index,:) + delta;
    % periodic boundary condition(PBC), verify whether particle is
    % out of the box
    % cdMoved(index,cdMoved(index,:)>L | cdMoved(index,:)< -L) = ...
    %     abs(abs(cdMoved(index,cdMoved(index,:)>L | cdMoved(index,:)<0)) - L);
    cdMoved(index,cdMoved(index,:)>L) = cdMoved(index,cdMoved(index,:)>L) - 2 * L;
    cdMoved(index,cdMoved(index,:)<-L) = cdMoved(index,cdMoved(index,:)<-L) + 2 * L;
    % position vetors with PBC
    r = cdMoved(index,:) - cdMoved;
    r = r - round(r/(2*L)).*(2*L);
    % calculate the norm of the position vectors
    norm_r = sqrt(sum(r.^2,2));
    norm_r(index) = [];
    % Figure Out Whether paritcle overlap
    if all(norm_r(:) > 2 * phys_cond.R)
        break;
    end
end
%initial position between index particle and the others
rInit = cdInit(index,:) - cdInit;
rInit = rInit - round(rInit/(2*L)).*(2*L);
norm_rInit = sqrt(sum(rInit.^2,2));
norm_rInit(index) = [];
%initial potential and Hamiltonian
HamInit = sum((1/(4*pi*phys_const.epsilon0*phys_const.epsilonr))*chargemulti./norm_rInit,"all");
% Hamiltonian after moving
Ham = sum((1/(4*pi*phys_const.epsilon0*phys_const.epsilonr))*chargemulti./norm_r,"all");
% difference between initial and after
delta_H = Ham - HamInit;
% Monte Carlo method below, determine the movement
% probability of moving
if exp(-(delta_H * 1e6)/(phys_const.kb * phys_cond.T))<= rand
    cdMoved = cdInit;
end
end