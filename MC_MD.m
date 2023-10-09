function MC_MD(Coordin_Char, L, T, N, k, k_b, R, C)
%%Initial Values
step = 200;
moved_times = 0;
moved_times_last = 0;
num_loops = 0;

while true
    index = ceil(N * rand);
    delta = step * normrnd(0, 1 , [1 , 3]);
    Coordin_Char_mvd = Coordin_Char;
    Coordin_Char_mvd(index , 1:3) = Coordin_Char(index , 1:3) + delta;

    %%whether it's out of the box
    for i = 1 : 3
        %         Coordin_Char_mvd(index , Coordin_Char_mvd(index,i)<0) = L(i) + Coordin_Char_mvd(index , Coordin_Char_mvd(index,i)<0);
        %         Coordin_Char_mvd(index , Coordin_Char_mvd(index,i)>L(i)) = Coordin_Char_mvd(index , Coordin_Char_mvd(index,i)>L(i)) - L(i);
        if Coordin_Char_mvd(index , i) < 0
            Coordin_Char_mvd(index , i) = L(i) + Coordin_Char_mvd(index , i);
        end
        if Coordin_Char_mvd(index , i) > L(i)
            Coordin_Char_mvd(index , i) = Coordin_Char_mvd(index , i) - L(i);
        end
    end

    posit_vect_mvd = Coordin_Char_mvd(: , 1:3) - Coordin_Char_mvd(index , 1:3);
    %posit_vect_mvd(all(a==0,2),:) = []; %remove lines with all elements are 0
    posit_vect_mvd = posit_vect_mvd - round(posit_vect_mvd ./ L) .* L;
    dist_mvd = sqrt(sum(posit_vect_mvd(: , :) .^ 2 , 2));

    %%whether there is something overlapped
    if (dist_mvd(:) < 2 * R) & (dist_mvd(:) ~= 0)
        continue;
    end

    %%Hamiltonian Calculation(initial)
    posit_vect = Coordin_Char(: , 1:3) - Coordin_Char(index , 1:3);
    posit_vect = posit_vect - round(posit_vect ./ L) .* L;
    dist = sqrt(sum(posit_vect(: , :) .^ 2 , 2));
    potential = ( k * Coordin_Char(: , 4) .* Coordin_Char(index , 4) ) ./ dist;
    potential(index) = []; %remove lines with all elements are 0
    Hamiltonian = sum(potential,'omitnan');

    %%Hamiltonian Calculation(moved)
    potential_mvd = ( k * Coordin_Char_mvd(: , 4) .* Coordin_Char_mvd(index , 4) ) ./ dist_mvd;
    potential_mvd(index) = [];
    Hamiltonian_mvd = sum(potential_mvd,'omitnan');

    delta_Ham = Hamiltonian_mvd - Hamiltonian;

    %%Boltzmann
    if delta_Ham <= 0
        Coordin_Char = Coordin_Char_mvd;
        moved_times = moved_times + 1;
    else%%probability of moving
        p = exp( -(delta_Ham * 1000) / (k_b * T) );
        if p <= rand
            Coordin_Char = Coordin_Char_mvd;
            moved_times = moved_times + 1;
        end
    end

    num_loops = num_loops + 1;
    if mod(num_loops,1e5) == 0
        if abs((moved_times - moved_times_last)/moved_times) < 0.005
            break;
        else
            moved_times_last = moved_times;
            moved_times = 0;
        end
    end
end

%% Write file
eval(['fileID = fopen(''coordinate_steady',num2str(C),'.txt''',',','''w''',');'])
fprintf(fileID,'%f\t %f\t %f\t %f\n',Coordin_Char');
fclose(fileID);
end