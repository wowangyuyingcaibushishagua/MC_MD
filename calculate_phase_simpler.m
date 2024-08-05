function [amplitude] = calculate_phase_simpler(phys_cond,cdMoved)
count = 1e3;
amplitude = 0;
% diff_phase = 0;
% center point of the box now is at the origin
% centerpoint = [0 0 0];
% cutoff distance
% rc = ceil(distances);
% rc = distances;
% generate photon several times 
for find_count = 1:count
    r = cdMoved;
    % Wave Vector k
    k = [1 0 0];
    % Decide eletric field E
    E = [0 0 1];
    % Determine whether the Point is at the Front of Wave
    % Remove Points BEHIND the Wave
    r(r(:,1)<=0,:) = [];
    % r = r - centerpoint;
    % distance of r(norm of r)
    norm_r = sqrt(sum(r.^2,2));
    % r(norm_r < (3 * rc),:) = [];
    % norm_r(norm_r < (3 * rc))=[];
    n1 = ceil(size(norm_r,1)*rand);
    r12 = r(n1,:) - r([1:(n1-1) (n1+1):end],:);
    d12 = sqrt(sum(r12 .* r12,2));
    % phasediff = 2*pi*(norm_r(n1)-norm_r(d12==max(d12))+abs(sum(r12(d12==max(d12),:).*k)))/phys_cond.wavelambda;
    E1 = cross(cross(r(n1,:),E),r(n1,:))/(norm_r(n1)^2) * exp(1i*2*pi*norm_r(n1)/phys_cond.wavelambda);
    E2 = cross(cross(r(d12==max(d12),:),E),r(d12==max(d12),:))/(norm_r(d12==max(d12))^2) * exp(1i*(2*pi*norm_r(d12==max(d12))+dot(r(d12==max(d12),:)-r(n1),k))/phys_cond.wavelambda);
    amplitude = amplitude + sum(real(E1 + E2).^2);
%     for n2 = 1:size(norm_r,1)
%         if ~(n2 == n1) && (norm_r(n2) > norm_r(n1))
%             r12 = (r(n2, :) - r(n1, :)) - round((r(n2, :) - r(n1, :)) / (2*L)) * (2*L);
%             d12 = sqrt(sum(r12 .* r12));
%             if d12 < (1.5 * rc)
%                 phasediff = 2*pi*(norm_r(n1)-norm_r(n2)+abs(sum(r12.*k)))/phys_cond.wavelambda;
%                 diff_phase = diff_phase + real(sqrt(sum(cross(r(n2,:),E).^2,"all"))/norm_r(n2)^2 + sqrt(sum(cross(r(n1,:),E).^2,"all"))/norm_r(n1)^2 * exp(1i*phasediff))^2;
%                 count = count + 1;
%             end
%         else
%             continue;
%         end
%     end
% end
% if count ~=0
%     diff_phase = diff_phase/count;
% else
%     diff_phase = 0;
% end
end
amplitude = amplitude / count;