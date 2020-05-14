%Run simulations
clear
clc

metabolites = {'GSH'};
[outA,outB,outC,outD] = sim_HERMES_Philips_Univ_SH(metabolites);

y_lim = [-0.1 0.20];
x_lim = [1.0 5];
figure(1),subplot(4,1,1),plot(outA.ppm,outA.specs,'r','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('A: 4.56ppm'),ylim(y_lim)
figure(1),subplot(4,1,2),plot(outB.ppm,outB.specs,'b','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('B: 1.90ppm'),ylim(y_lim)
figure(1),subplot(4,1,3),plot(outC.ppm,outC.specs,'g','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('C: 4.56ppm 1.90ppm'),ylim(y_lim)
figure(1),subplot(4,1,4),plot(outD.ppm,outD.specs,'k','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('D: non-editing'),ylim(y_lim)

% Figure for GABA On 
%figure(2), plot(outA.ppm,outB.specs+outC.specs-outA.specs-outD.specs,'b',outA.ppm,outA.specs-outD.specs,'r',outA.ppm,outB.specs-outC.specs,'k'),set(gca,'xdir','reverse'),xlim([1 4]), xlabel('ppm'),legend('GABA On','A-D','B-C')

% Figure for GSH On
figure(3), plot(outA.ppm,outA.specs+outC.specs-outB.specs-outD.specs,'b',outA.ppm,outB.specs+outC.specs-outA.specs-outD.specs,'r'),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),legend('GSH On','GABA Off')