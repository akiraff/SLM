%% compare two SLM pattern
T13 = readmatrix('SLM_phase_ring_screen_13Bit.csv');
%%
T14 = readmatrix('SLM_phase_screen_ring_14Bit.csv');
%%
figure(100);clf;
subplot(1,2,1);
imagesc(T13);
axis equal;
colorbar
title('13 Bit')
subplot(1,2,2);
imagesc(T14);
axis equal;
colorbar
title('14 Bit');
%%
figure(101);clf;
Tdiff = T13-T14;
imagesc(Tdiff);
axis equal
colorbar;
%%
P13 = T13(500,:);
P14 = T14(500,:);
figure(200);clf;
plot(1:length(P13), P13, 'ro', 1:length(P14), P14, 'bo');
legend('13 Bit','14 Bit');
title('Row 500');