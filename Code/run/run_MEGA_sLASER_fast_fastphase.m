%Run HERCUsLASER simulations

clc
close all

tic

addpath /Users/steve/Documents/My_Studies/Hercules_2_Study/Code/
addpath /Users/steve/Documents/My_Studies/Hercules_2_Study/Code/Simulation/Functions/sim_sLASER_based/

% run GOIA_build_sweep_v4 to generate RF_GOIA.mat and Sweep
metabolites = {'GABA'};
[outA,outB] = sim_MEGA_sLASER_fast_fastphase_SH(metabolites);
toc
%time_elapsed=toc;
% outSum.ppm = outA.ppm;
% outSum.specs = (outA.specs+outB.specs+outC.specs+outD.specs)/4;
% Lac_DIFF = Lac_fast_fast_outB.specs+Lac_fast_fast_outD.specs-Lac_fast_fast_outA.specs-Lac_fast_fast_outC.specs;
% Lac_resid = Lac_fast_fast_outA.specs+Lac_fast_fast_outB.specs-Lac_fast_fast_outC.specs-Lac_fast_fast_outD.specs;
% figure(3), plot(Lac_fast_fast_outA.ppm,Lac_DIFF,'b',Lac_fast_fast_outA.ppm,Lac_resid,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('Lac Diff','Lac Residual')
% 
% 
x_lim = [1 5];
%y_lim = [-0.2 0.4];
figure(1),subplot(4,1,1),plot(outA.ppm,outA.specs,'r','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('A: 7.5ppm')
figure(1),subplot(4,1,2),plot(outB.ppm,outB.specs,'b','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('B: 1.9ppm')

% plot GABA
figure(2), plot(outA.ppm,outB.specs-outA.specs,'k'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('GABAGlx')
% GABA_DIFF = outC.specs+outD.specs-outA.specs-outB.specs;
% GABA_resid = outA.specs+outC.specs-outB.specs-outD.specs;
% figure(3), plot(outA.ppm,GABA_DIFF,'b',outA.ppm,GABA_resid,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('GABA Diff','GABA Residual')
% 
% % plot GSH
% % figure(4), plot(outA.ppm,outA.specs+outC.specs-outB.specs-outD.specs),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('GSH')
% % 
% % GSH_DIFF = outA.specs+outC.specs-outB.specs-outD.specs;
% % GSH_resid = outC.specs+outD.specs-outA.specs-outB.specs;
% % GSH_DIFF_new = new_Lac.outA.specs+new_Lac.outC.specs-new_Lac.outB.specs-new_Lac.outD.specs;
% % GSH_resid_new = new_Lac.outC.specs+new_Lac.outD.specs-new_Lac.outA.specs-new_Lac.outB.specs;
% % 
% % figure(5), plot(outA.ppm,GSH_DIFF/4,'b',outA.ppm,GSH_resid/4,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('old','new')
% % figure(6), plot(outA.ppm,GSH_resid/4,'b',outA.ppm,GSH_resid_new/4,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('old','new')
% % 
% % figure(7), plot(outA.ppm,GSH_DIFF_new/4,'b',outA.ppm,GSH_resid_new/4,'r'),set(gca,'xdir','reverse'),xlim(x_lim),ylim(y_lim), xlabel('ppm'),legend('old','new')
