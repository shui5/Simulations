% sim_press_shaped_fastphase.m
% Adopted from sim_press_shaped.m from Jamie Near
% Steve C.N. Hui, Muhammad G. Saleh and Richard A.E. Edden. 2020
% The Johns Hopkins University School of Medicine
% 
% USAGE:
% out = sim_press_shaped_fastphase(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,dy,Gx,Gy,phCyc1,phCyc2,flipAngle,centreFreq)
% 
% DESCRIPTION:
%

function out = sim_PRESS_shaped_fastphase(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,dy,Gx,Gy,phCyc1,phCyc2,flipAngle,centreFreq)

if nargin<17
    centreFreq=3.0;%2.3/4.65;
    if nargin<16
        flipAngle=180;
    end
end
    
if tau1<tp/1000
    error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
end
if tau2<tp/1000
    error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
end

%Set water to centre
%centreFreq=2.3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from tau1 and tau2;
delays=zeros(2);
delays(1)=tau1-tp;
delays(2)=tau2-tp;
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end

d_temp=cell(length(phCyc1));

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                    %EXCITE
d=sim_evolve(d,H,delays(1)/2000);                         %Evolve by delays(1)/2
d_in=d;

for kk=1:length(phCyc1)
d_temp{kk}=sim_shapedRF(d_in,H,RF,tp,flipAngle,phCyc1(kk),dx,Gx);          %1st shaped 180 degree refocusing pulse
% average as we go
end

receiver_phase = exp(1i*phCyc1/180*pi*2);
%Copy from FID-A , calculate the average density matrix (Doing this inside a separate for 
%loop because I couldn't figure out how to do this inside the parfor loop): 
for kk=1:length(phCyc1)
   if kk==1
       d_out1=d_temp{kk};
   else%for RP1=1:length(refPhCyc1)
        d_out1=sim_dAdd(d_out1,d_temp{kk},receiver_phase(kk));
   end
    %end
end

d=sim_evolve(d_out1,H,(delays(1)+delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
%end of first phase cycle loop scnh
clear d_temp;

d_in=d;
for kk=1:length(phCyc2)
d_temp{kk}=sim_shapedRF(d_in,H,RF,tp,flipAngle,phCyc2(kk),dy,Gy);          %1st shaped 180 degree refocusing pulse
end

receiver_phase = exp(1i*phCyc2/180*pi*2);
%Copy from FID-A , calculate the average density matrix (Doing this inside a separate for 
%loop because I couldn't figure out how to do this inside the parfor loop): 
for kk=1:length(phCyc2)
   if kk==1
       d_out2=d_temp{kk};
   else%for RP1=1:length(refPhCyc1)
        d_out2=sim_dAdd(d_out2,d_temp{kk},receiver_phase(kk));
   end
    %end
end

%d=sim_shapedRF(d,H,RF,tp,flipAngle,90+phCyc2,dy,Gy);          %2nd shaped 180 degree refocusing pulse
%ned of second phase cycle loop scnh

d=sim_evolve(d_out2,H,delays(2)/2000);                            %Evolve by delays(2)/2
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='press';
out.te=tau1+tau2;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isISIS=0;








 
