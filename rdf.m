function [r,g,r_at_maxg] = rdf(Coordin, L, pbc, Ng, N)
L_times_pbc = L .* pbc; % deal with boundary conditions
rho = N / prod(L);      % global particle density
rc = sqrt(L(1)^2+L(2)^2+L(3)^2)/2;  % the maximum radius
dr = rc / Ng;           % bin size
%% accumulate
g = zeros(Ng, 1);
for n1 = 1 : (N - 1)                                % sum over the atoms
    for n2 = (n1 + 1) : N                           % skipping half of the pairs
        r12 = Coordin(n2, :) - Coordin(n1, :);                  % position difference vector
        r12 = r12 - round(r12 ./L ) .* L_times_pbc; % minimum image convention
        d12 = sqrt(sum(r12 .* r12));                % distance
        if d12 < rc                                 % there is a cutoff
            index = ceil(d12 ./ dr);                % bin index
            g(index) = g(index) + 1;                % accumulate
        end
    end
end
%% normalize
for n = 1 : Ng
    g(n) = g(n) / N * 2;           % 2 because half of the pairs have been skipped
    dV = 4 * pi * (dr * n)^2 * dr; % volume of a spherical shell
    g(n) = g(n) / dV;              % now g is the local density
    g(n) = g(n) / rho;             % now g is the RDF
end
%%
[~,index] = max(g);
r_at_maxg = index * dr;%r corresponding max g
r = transpose((1 : Ng) * dr);
g = smooth(smooth(smooth(g)));
end