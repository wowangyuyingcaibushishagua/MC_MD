function [cdInit,charge,numParticle] = Init(C,phys_cond,phys_const)
numParticle = floor(C * phys_const.Na / 1e13 * phys_cond.L^3);
cdInit = phys_cond.L * rand(numParticle,3) /10;
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
charge = phys_const.e * poissrnd(phys_cond.averagee,numParticle,1);
end