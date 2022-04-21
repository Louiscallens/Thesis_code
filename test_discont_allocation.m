for i = 1:20
N = 20;
uniform = linspace(0, 1, 20);

nb_disconts = 4;
disconts = sort(rand(1,4));

figure; clf;
subplot(211); hold off;
xline(uniform, 'b'); hold on;
xline(disconts, 'r');

modified = init_intervals_discontinuities(disconts, N, 1);

subplot(212); hold off;
xline(modified, 'b'); hold on;
xline(disconts, 'r');

disp(length(uniform));
disp(length(modified));
pause();
end