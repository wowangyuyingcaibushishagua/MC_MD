function [diff_phase] = calculate_phase(phys_cond,cdMoved,L)
count = 0;diff_phase = 0;
% center point of the box
centerpoint = [0 0 0];
% cutoff distance
% rc = ceil(distances);
rc = 1000;
% generate photon several times 
for gen_time = 1:1e4
    r = cdMoved;
    % generate wave vector randomly in spherical coordinates
    %Function 'rand()' generates random number in range of 0 to 1 uniformly
    theta = pi * rand;
    phi = 2 * pi * rand;
    if sin(theta) * sin(phi) == 0 && sin(theta) * cos(phi) == 0 && cos(theta) ==0
        continue;
    end
    % Wave Vector k
    k = [sin(theta) * sin(phi) sin(theta) * cos(phi) cos(theta)];
    % Decide eletric field E
    if cos(theta) ~= 0
        Ex = rand;Ey = rand;
        E = [Ex Ey -(sin(theta) * sin(phi) * Ex + sin(theta) * cos(phi) * Ey)/cos(theta)];
    else
        Ez = rand;
        if sin(theta) * cos(phi) ~= 0
            Ex = rand;
            E = [Ex -(sin(theta) * sin(phi) * Ex)/(sin(theta) * cos(phi))];
        else
            E = [0 rand Ez];
        end
    end
    % normalization of E
    E = E./sqrt(sum(E.^2,"all"));
    % Determine whether the Point is at the Front of Wave
    % Angle between Position Vector and Wave Vector
    rdotn = sum((cdMoved - centerpoint) .* k,2);
    % Remove Points BEHIND the Wave
    r(rdotn(:)<0,:) = [];
    % r = (r - centerpoint) - round((r - centerpoint)/L).*L;
    r = (r - centerpoint);
    % distance of r(norm of r)
    norm_r = sqrt(sum(r.^2,2));
    % r(norm_r < (4 * rc) | norm_r > (6 * rc),:) = [];
    % norm_r(norm_r < (4 * rc) | norm_r > (6 * rc))=[];
    % r(norm_r(:) < (1.8 * rc) | norm_r(:) > (3 * rc),:) = [];
    % norm_r(norm_r(:) < (1.8 * rc) | norm_r(:) > (3 * rc))=[];
    r(norm_r < (3 * rc),:) = [];
    norm_r(norm_r < (3 * rc))=[];
    % % sum over the atoms
    % for n1 = 1 : size(r,1)
    %     % skipping half of the pairs
    %     for n2 = (n1 + 1) : size(r,1)
    %         % position vector between particles
    %         r12 = r(n2, :) - r(n1, :);
    %         % periodic boundary condition
    %         r12 = r12 - round(r12 / L) * L;
    %         % distance r_ij
    %         d12 = sqrt(sum(r12 .* r12));
    %         if d12 < (3 *rc)
    %             phasediff = 2*pi*abs(norm_r(n1)-norm_r(n2)+abs(sum(r12.*n)))/phys_cond.wavelambda;
    %             diff_phase = diff_phase + (3/(16*pi) * (abs(sum(r(n1,:).*n)/norm_r(n1))^2) * real(1/norm(n2) + 1/norm(n1) * exp(1i*phasediff)))^2;
    %             count = count + 1;
    %         end
    %     end
    % end
    n1 = ceil(size(norm_r,1)*rand);
    for n2 = 1:size(norm_r,1)
        if ~(n2 == n1) && (norm_r(n2) > norm_r(n1))
            % r12 = (r(n2, :) - r(n1, :)) - round((r(n2, :) - r(n1, :)) / L) * L;
            r12 = (r(n2, :) - r(n1, :)) - round((r(n2, :) - r(n1, :)) / (2*L)) * (2*L);
            d12 = sqrt(sum(r12 .* r12));
            % if d12 < (1.2 * rc) || d12 > (0.8 * rc)
            if d12 < (3 *rc)
                phasediff = 2*pi*(norm_r(n1)-norm_r(n2)+abs(sum(r12.*k)))/phys_cond.wavelambda;
                diff_phase = diff_phase + real(sqrt(sum(cross(r(n2,:),E).^2,"all"))/norm_r(n2)^2 + sqrt(sum(cross(r(n1,:),E).^2,"all"))/norm_r(n1)^2 * exp(1i*phasediff))^2;
                count = count + 1;
            end
        else
            continue;
        end
    end
end
if count ~=0
    diff_phase = diff_phase/count;
else
    diff_phase = 0;
end
end