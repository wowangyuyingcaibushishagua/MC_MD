function [cdInit,charge,L] = Init_C(C,phys_cond,phys_const)
% numParticle = floor(C * phys_const.Na / 1e13 * phys_cond.L^3);
L = 0.5 * (1e13 * phys_cond.numParticles / (phys_const.Na * C)) ^ (1/3);
cdInit = L * rand(phys_cond.numParticles,3) * (-1)^(round(rand));
% cdInit = zeros(numParticle,3);
% count = 1;
% for i = 1:ceil(numParticle ^ (1 / 3))
%     for j = 1:ceil(numParticle ^ (1 / 3))
%         for k = 1:ceil(numParticle ^ (1 / 3))
%             cdInit(count,1) = 500 + (i-1)*500;
%             cdInit(count,2) = 500 + (j-1)*500;
%             cdInit(count,3) = 500 + (k-1)*500;
%             count = count+1;
%             if count == numParticle+1
%                 break
%             end
%         end
%         if count == numParticle+1
%             break
%         end
%     end
%     if count == numParticle+1
%         break
%     end
% end
charge = phys_const.e * poissrnd(phys_cond.averagee,phys_cond.numParticles,1);
end