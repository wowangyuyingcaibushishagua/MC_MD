function [g] = rdf(cdMoved, L, Ng, numParticle)
% cutoff, should change with the concentration or D
rc = 20000;
% global particle density
rho = numParticle / ((2*L)^3);
% bin size
dr = rc / Ng;
% accumulation
g = zeros(Ng, 1);
for i = 1 : (numParticle - 1)% sum over the atoms
    for j = (i + 1) : numParticle% skipping half of the pairs
        % position difference vector
        r12 = cdMoved(j, :) - cdMoved(i, :);
        % minimum image convention
        r12 = r12 - round(r12 ./ (2*L)) .* (2*L);
        % distance
        d12 = sqrt(sum(r12 .* r12));
        if d12 < rc% there is a cutoff
            % bin index
            index = ceil(d12 / dr);
            % accumulation
            g(index) = g(index) + 1;
        end
    end
end
% normalization
for i = 1 : Ng
    % 2 because half of the pairs have been skipped
    g(i) = g(i) / numParticle * 2;
    % volume of a spherical shell
    dV = 4 * pi * (dr * i)^2 * dr;
    % now g is the local density
    g(i) = g(i) / dV;
    % now g is the RDF
    g(i) = g(i) / rho;
end