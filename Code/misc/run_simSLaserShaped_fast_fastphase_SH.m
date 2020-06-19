% run_simSLaserShaped_fast_fastphase.m

% Fast sLASER with fast phase cycling version by 
% Steve C.N. Hui, Muhammad G Saleh and Richard A.E. Edden
% (The Johns Hopkins University School of Medicine, 2020).

% Adopted from run_simSLaserShaped_fast.m
% Dana Goerzen and Jamie Near (McGill University, 2019).
%
% USAGE:
%
% DESCRIPTION:
% Simulates CMRR semi-LASER experiments designed by Oz et al. (2018) using
% shaped adiabatic refocusing pulses along XX YY gradients
% The excitation is simulated as an instantaneous rotation,
% and the adiabatic refocusing pulses are simulated as a shaped rotation.
% This simulation also employs a 4-step phase cycling scheme as follows:

% signal = ([0 0, 0 0] - [0 0, 0 90]) - ([0 90, 0 0] + [0 90, 0 90]);

%
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the sLASER experiment, you have to run this
% simulation many times at various points in space (x,y), and then add
% together the resulting spectra and scale down by the number of simulations.
%
% Feb 2020 - Jamie Near:  This code now accepts gradient modulated pulses.  
%

%tic;

% INPUTS:
n=8192;%4096; %= number of points in fid/spectrum
sw=4000;%5000; %= desired spectral width in [Hz]
Bfield=3;%7; %= main magnetic field strength in [T]
lw=1;%2; %= linewidth in [Hz]
load spinSystems.mat; %= spin system definition structure
sys=sysH2O;
%rfPulse=io_loadRFwaveform('sampleAFPpulse_HS2_R15.RF','inv'); % adiabatic RF pulse shaped waveform
load RF_GOIA_20200506_100pts.mat; % adiabatic RF pulse shaped waveform
rfPulse=Sweep2;

refTp=4.5; %4.5 (for GOIA);%9.395 (for sampleAFPpulse_HS2_R15);%3.5; %= RF pulse duration in [ms]
flipAngle=180; %= flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
centreFreq=3.0;%4.65;%2.3; %= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
thkX=3;%2; %slice thickness of x refocusing pulse [cm]
thkY=3;%2; %slice thickness of y refocusing pulse [cm]
fovX=4.5;%3; %size of the full simulation Field of View in the x-direction [cm]
fovY=4.5;%3; %size of the full simulation Field of View in the y-direction [cm]
nX=41;%16; %Number of grid points to simulate in the x-direction
nY=41;%16; %Number of grid points to simulate in the y-direction
te=30;%135;         %sLASER total echo time [ms]
% ph1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse
% ph2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse
% ph3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse
% ph4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse

ph1=[0 90 180 270]; %phase cycling scheme of first refocusing pulse
ph2=[0 90 180 270]; %phase cycling scheme of second refocusing pulse
ph3=[0 90 180 270]; %phase cycling scheme of third refocusing pulse
ph4=[0 90 180 270]; %phase cycling scheme of fourth refocusing pulse

% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS
%             sequence.

%set up spatial grid
% x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
% y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]
if nX==1 
x=0;
else
    x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
end

if nY==1
    y=0;
else
    y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]
end

gamma=42577000; %gyromagnetic ratio

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
%rfPulse=rf_resample(rfPulse,100);

%sys=sysRef0ppm
if ~rfPulse.isGM
    %Non-gradient modulated pulse - Calculating the x and y gradient 
    %strengths for the desired slice thickness
    Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
else
    %Gradient modulated pulse
    %1.  Calculating the unitless scaling factor for the GM waveform.
    Gx=(rfPulse.tthk/(refTp/1000))/thkX; % added a pair of bracket, scnh
    Gy=(rfPulse.tthk/(refTp/1000))/thkY; % added a pair of bracket, scnh
end

%Initialize structures:
% out_posxy_rpc=cell(length(x),length(y),length(ph1));

%out_posx_rpc =cell(length(x),length(ph1)); % can we remove this? scnh
out_posx =cell(length(x)); % added. scnh
%d=cell(length((ph1)));
d=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
parfor X=1:length(x)
%for X=1:length(x)
    %for m=1:length(ph1)
    disp(['Executing X-position ' num2str(X)]);
    %    ...
    %        '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
    %out_posx_rpc{X}{m}=sim_sLASER_shaped_Ref1(Bfield,sys,te,rfPulse,refTp,x(X),Gx,ph1(m),ph2(m),flipAngle,centreFreq);
    out_posx{X}=sim_sLASER_shaped_Ref1_fastphase(Bfield,sys,te,rfPulse,refTp,x(X),Gx,ph1,ph2,flipAngle,centreFreq);
end

for X=1:length(x)
    %for m=1:length(ph1)
    d=sim_dAdd(d,out_posx{X});
    %end
end

% %Initialize structures:
out_posy =cell(length(y));
out=struct([]);

%Now loop through y direction (second refoc pulse only);
for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
%parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
    %for m=1:length(ph1)
    disp(['Executing Y-position ' num2str(Y)]);
    %        ...
    %            '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
    out_posy{Y}=sim_sLASER_shaped_Ref2_fastphase(d,n,sw,Bfield,lw,sys,te,rfPulse,refTp,y(Y),Gy,ph3,ph4,flipAngle,centreFreq);
    %end
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
    %for m=1:length(ph1)
        out=op_addScans(out,out_posy{Y}); % should the function be sim_dAdd or op_addScans? scnh
    %end
end


%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY*(length(ph1)*length(ph2)*length(ph3)*length(ph4)));
out=op_ampScale(out,1/numSims);
%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);
%Plot output
figure(1),plot(out.ppm,out.specs*exp(1i*-0*pi/180)),xlim([0 4.5])

%time_sim=toc

%%%%%NESTED FUNCTIONS BELOW%%%%%%%%%

%% Simulate in X-direction only
function d = sim_sLASER_shaped_Ref1_fastphase(Bfield,sys,te,RF,tp,dx,Gx,ph1,ph2,flipAngle,centreFreq)

%if nargin<21
    %centreFreq=2.3;%4.65;%2.3;
%     if nargin<20
%         flipAngle=180;
%     end
%end

%Check if this is a gradient modulated pulse.  If so, set Gx equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gx and then set Gx equal to zero:
    RF=rf_scaleGrad(RF,Gx);
    Gx=0;
end

if (te/4)<(tp/1000)
    error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
end

%initialize evolution times
tau1=(te/4-tp)/2;
tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d_in=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
% loop over ph1
for kk=1:length(ph1)
    d_out1{kk}=sim_shapedRF(d_in,H,RF,tp,flipAngle,ph1(kk),dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
    
end
% end loop over ph1 and average over ph1
% for  kk=1:length(ph1)
%     d_out1{kk}
% end
%out_1=struct([]);

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
%for Y=1:length(y)
%for RP1=1:length(refPhCyc1)
for RP1=1:length(ph1)
    receiver_phase=round((exp(1i*ph1(RP1)/180*pi*2)));
    %out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
    %out_1=sim_dAdd(out_1,d_out1{RP1},receiver_phase);
    if RP1==1
        out_1=d_out1{RP1};
    else
        out_1=sim_dAdd(out_1,d_out1{RP1},receiver_phase);
    end
end
%end
%end
%clear d_in;
d_in=sim_evolve(out_1,H,tau2/1000);                            %Evolve by tau2

for kk=1:length(ph2)
    d_out2{kk}=sim_shapedRF(d_in,H,RF,tp,flipAngle,ph2(kk),dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
    
end
% end loop over ph1 and average over ph1
% for  kk=1:length(ph1)
%     d_out1{kk}
% end
%out_2=struct([]);

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
%for Y=1:length(y)
%for RP1=1:length(refPhCyc1)
for RP2=1:length(ph2)
    receiver_phase=round((exp(1i*ph2(RP2)/180*pi*2)));
    %out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
    %out_2=sim_dAdd(out_2,d_out2{RP2},receiver_phase);
    
    if RP2==1
        out_2=d_out2{RP2};
    else
        out_2=sim_dAdd(out_2,d_out2{RP2},receiver_phase);
    end
    
end
        

d=sim_evolve(out_2,H,tau2/1000);       %should d be out_2? scnh                     %Evolve by tau2

end

%% Simulate in Y-direction only
function out = sim_sLASER_shaped_Ref2_fastphase(d,n,sw,Bfield,linewidth,sys,te,RF,tp,dy,Gy,ph3,ph4,flipAngle,centreFreq)

% if nargin<21
%     centreFreq=2.3;%4.65;%2.3;
%     if nargin<20
         flipAngle=180;
%     end
% end

%Check if this is a gradient modulated pulse.  If so, set Gy equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gy and then set Gy equal to zero:
    RF=rf_scaleGrad(RF,Gy);
    Gy=0;
end

if (te/4)<(tp/1000)
    error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
end

%initialize evolution times
tau1=(te/4-tp)/2;
tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%d_out3=cell([]);
%BEGIN sLASER PULSE SEQUENCE************
% loop over ph1 
for kk=1:length(ph3)
    %d_out1{kk}=sim_shapedRF(d_in,H,RF,tp,flipAngle,ph1(kk),dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
    d_out3{kk}=sim_shapedRF(d,H,RF,tp,flipAngle,ph3(kk),dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
    
end

%out_3=struct([]);

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
%for Y=1:length(y) 
    %for RP1=1:length(refPhCyc1)
        for RP3=1:length(ph3)
            receiver_phase=round((exp(1i*ph3(RP3)/180*pi*2)));
            %out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
             %out_3=sim_dAdd(out_3,d_out3{RP3},receiver_phase);
             
             if RP3==1
                 out_3=d_out3{RP3};
             else
                 out_3=sim_dAdd(out_3,d_out3{RP3},receiver_phase);
             end
             
        end
    %end
%end
d=sim_evolve(out_3,H,tau2/1000);                            %Evolve by tau2

for kk=1:length(ph4) 
d_out4{kk}=sim_shapedRF(d,H,RF,tp,flipAngle,ph4(kk),dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
end

%out_4=struct([]);

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:

for RP4=1:length(ph4)
    receiver_phase=round((exp(1i*ph4(RP4)/180*pi*2)));
    %out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
    %out_4=sim_dAdd(out_4,d_out4{RP4},receiver_phase);
    if RP4==1
        out_4=d_out4{RP4};
    else
        out_4=sim_dAdd(out_4,d_out4{RP4},receiver_phase);
    end
end

d=sim_evolve(out_4,H,tau1/1000);                            %Evolve by tau1

[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
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
end
