function initialization(N, C)
r = zeros(N,3);
count = 1;
for i = 1:ceil(N ^ (1 / 3))
    for j = 1:ceil(N ^ (1 / 3))
        for k = 1:ceil(N ^ (1 / 3))
            r(count,1) = 500 + (i-1)*500;
            r(count,2) = 500 + (j-1)*500;
            r(count,3) = 500 + (k-1)*500;
            count = count+1;
            if count == N+1
                break
            end
        end
        if count == N+1
            break
        end
    end
    if count == N+1
        break
    end
end
%% print coordinates to file
eval(['fid= fopen(''coordinate',num2str(C),'.txt''',',','''w''',');'])
for i = 1 : N
    e = 1.6 * poissrnd(9.6);
    fprintf(fid, '%f\t', r(i,1));
    fprintf(fid, '%f\t', r(i,2));
    fprintf(fid, '%f\t', r(i,3));
    fprintf(fid, '%f\r', e);
end
fclose(fid);
end