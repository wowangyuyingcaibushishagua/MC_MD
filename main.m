clear all;
%% Physical Constans
epsilon_r = 80;
k = 8.9876 / epsilon_r;%Column Constant
k_b = 1.3806;%Boltzmann constant
%% Experimental Conditions
T = 300;
L = [21544 21544 21544];%Box Size
R = 7.5;%Radius of AuNPs
pbc = [1, 1, 1];%periodic boundary condition
%C = 65;%unit is in pM for pmol/L
r_at_maxg = zeros();
for C = 6 : 1 : 65
    %%
    N = num_of_particles(C);%from Molarity get number of particles in Box
    %%
    initialization(N, C);
    %%
    eval(['Coordin_Char = load(''coordinate',num2str(C),'.txt''',');'])
    %% distribution initial state
    % figure(1)
    % for i=1:N
    %     Plot_Spheres(Coordin_Char(i,1),Coordin_Char(i,2),Coordin_Char(i,3),R);
    % %     hold on
    % %     axis equal
    %     grid on
    % end
    %%
    MC_MD(Coordin_Char, L, T, N, k, k_b, R, C);
    %%
    eval(['Coordin_Char_steady = load(''coordinate_steady',num2str(C),'.txt''',');'])
    %% distribution steady state
    % figure(2)
    % for i=1:N
    %     Plot_Spheres(Coordin_Char_steady(i,1),Coordin_Char_steady(i,2),Coordin_Char_steady(i,3),R*10);
    % %     hold on
    % %     axis equal
    %     grid on
    % end
    %% number of bins (number of data points in the figure below)
    Ng = 900;
    %%
    [r,g,r_at_maxg(C)] = rdf(Coordin_Char_steady(:, 1 : 3), L, pbc, Ng, N);
    %% plot rdf
    % figure(3);
    % plot(r, g, '-');
    % xlim([0, 12000]);
    % ylim([0, 1.5]);
    % xlabel('r /nm', 'fontsize', 15);
    % ylabel('g(r)', 'fontsize', 15);
    % set(gca, 'fontsize', 15);
    % grid on
end
C = 1:1:65;
figure;
plot(C,r_at_maxg,'o-');
xlabel('C/pM', 'fontsize', 15);
ylabel('Corresponding Distance/nm', 'fontsize', 15);
xlim([6, 65]);
ylim([2000, 6000]);