%Run simulations"
metabs = {'GSH'};%,'GSH'};
%freq_ppms={[1.9 7.5]};%,[4.56 7.5]};
[outON,outOFF] = sim_MEGA_PRESS_Philips_Univ_SH(metabs);%,freq_ppms);
figure(1),plot(outON.ppm,outON.specs,'r',outON.ppm,outOFF.specs,'b',outON.ppm,(outON.specs-outOFF.specs)-0.2*ones(size(outOFF.specs)),'b'),set(gca,'xdir','reverse'),xlim([0 5]), xlabel('ppm')
