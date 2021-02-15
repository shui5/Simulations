% sim_MEGAsLASER_fast_SH.m

% Written by Steve CN Hui and Muhammad G Saleh, Johns Hopkins University, 2019

% Dana Goerzen, Jamie Near, McGill University, 2019
%
% DESCRIPTION:
% This script simulates experiments using MEGA sLASER localization with fully shaped editing
% and refocusing pulses. Phase cycling of both the editing and refocusing
% pulses is performed. Furthermore, simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal. Summation across phase cycles and spatial positions is
% performed. As a result of the phase cycling and spatially resolved simulations,
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations. To enable the use
% of the MATLAB parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes. If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
%
% Parameters description:
%
% refocWaveform     = name of refocusing pulse waveform.
% editWaveform1      = name of editing pulse waveform.
% editOnFreq        = freqeucny of edit on pulse[ppm]
% editOffFreq       = frequency of edit off pulse[ppm]
% refTp             = duration of refocusing pulses[ms]
% editTp            = duration of editing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% thkZ              = slice thickness of z excitation pulse [cm]
% x                 = vector of x positions to simulate [cm]
% y                 = vector of y positions to simulate [cm]
% z                 = vector of z positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]

function [outA,outB] = sim_MEGA_sLASER_fast_fastphase_SH_UQL(metabolites)

for kk=1:numel(metabolites)
    metabolite = metabolites{kk};
    % ********************PARAMETERS**********************************
    
    % Spin system to simulate
    spinSys     = metabolite;
    out_name    = ['Siemens_7T_MEGA_sLASER_101pts_' spinSys '.mat'];
    
    % General properties of the simulations
    Npts    = 8192;     % number of spectral points
    sw      = 4000;     % spectral width [Hz]
    Bfield  = 6.9809;        % Philips magnetic field strength [Tesla]
    lw      = 1;        % linewidth of the output spectrum [Hz]
    gamma   = 42577000; % gyromagnetic ratio
    
    % Define the pulse waveforms here
    %exciteWaveform      = 'sg100_200pts.pta';            % name of excitation pulse waveform.
    refWaveform          = 'HSN_waveform'; % adiabatic RF pulse shaped waveform
    editWaveform         = 'sl_univ_pulse.pta';   % name of 1st single editing pulse waveform. %ON/OFF
    
    % Define frequency parameters for editing targets
    flipAngle   = 180;
    centreFreq  = 3.0;          % Center frequency of MR spectrum [ppm]
    editOnFreq1 = 1.9;          % 4.56 for GSH, 1.9 for GABA     % Center frequency of 1st sLASER_PRESS experiment;
    editOnFreq2 = 7.5;          % Center frequency of 2nd sLASER_PRESS experiment;
    
    % Define pulse durations and flip angles specific for every vendor
    refTp       = 6;              % duration of refocusing pulses [ms], set to zero for hard pulse
    editTp1     = 10.24;          % duration of 1st editing pulse [ms]
    TE          = 74;               % Echo time [ms]
    
    % Define phase cycling pattern - stays the same, there are two editing
    % pulses in MEGA-PRESS with the same frequency, different time
    editPhCyc1=[0 90 180 270]; %phase cycling steps for 1st editing pulse [degrees]
    editPhCyc2=[0 90 180 270]; %phase cycling steps for 2nd editing pulse [degrees] From Shaped edit
    ph1=[0 90 180 270];  %phase cycling scheme of first refocusing pulse
    ph2=[0 90 180 270]; %phase cycling scheme of second refocusing pulse
    ph3=[0 90 180 270]; %phase cycling scheme of third refocusing pulse
    ph4=[0 90 180 270]; %phase cycling scheme of fourth refocusing pulse
    
    % Define spatial resolution of simulation grid
    fovX    = 4.5;  % size of the full simulation Field of View in the x-direction [cm]
    fovY    = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
    fovZ    = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
    thkX    = 3.0;    % slice thickness of x refocusing pulse [cm]
    thkY    = 3.0;    % slice thickness of y refocusing pulse [cm]
    thkZ    = 3.0;    % slice thickness of z excitation pulse [cm]
    nX      = 101;   % number of spatial points to simulate in x direction
    nY      = 101;   % number of spatial points to simulate in y direction
    nZ      = 101;   % number of spatial points to simulate in z direction
    if nX>1
        x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
    else
        x=0;
    end
    if nY>1
        y=linspace(-fovY/2,fovY/2,nY); %X positions to simulate [cm]
    else
        y=0;
    end
    if nZ>1
        z=linspace(-fovZ/2,fovZ/2,nZ); %X positions to simulate [cm]
    else
        z=0;
    end
    
    % set up exact timing - based on 'normal' pulse set on Philips 3T - SH
    %taus = [5.0688, abs(5.0688 - 24.3619), (38.3882-24.3619), (43.0007-38.3882), (49.6813-43.0007), (64.3619-49.6813), (80.0-64.3619)];
%    taus = [4.7523, (23.9861-4.7523), (38.1735-23.9861), (42.8285-38.1735), (49.4073-42.8285), (63.9861-49.4073), (80.0-63.9861)];
 taus = [(5.250-3.53)+((2.320+10.240)/2)+2.37,...
     (17.210-15.490)+((6.000+10.240)/2),...
     (26.650-23.210)+((6.000+6.000)/2),...
     (36.090-32.650)+((6.000+6.000)/2),...
     (43.810-42.090)+((6.000+10.240)/2),...
     (63.650-54.050)+((6.000+10.240)/2),...
     (74.0-69.650)+(6.000/2)];
    % ********************SET UP SIMULATION**********************************
    
    % Load RF waveforms for excitation and refocusing pulses
    
    %excRF           = io_loadRFwaveform(exciteWaveform,'exc',0);
    switch refWaveform
        case 'HSN_waveform'
            load hsecant.mat;  %Fill this up
            refRF      = hsecant;
    end
    
    if ~refRF.isGM
        %Non-gradient modulated pulse - Calculating the x and y gradient
        %strengths for the desired slice thickness
        Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
        Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
    else
        %Gradient modulated pulse
        %1.  Calculating the unitless scaling factor for the GM waveform.
        Gx=(refRF.tthk/(refTp/1000))/thkX;
        Gy=(refRF.tthk/(refTp/1000))/thkY;
    end
    
    editRF_A          = io_loadRFwaveform(editWaveform,'inv',0); %4.58
    %editRF_A.tw1      = 0.93;   % acquired from check_tw1.m
    % Construct the editing pulses from the waveforms and defined
    % frequencies
    editRFonA       =rf_freqshift(editRF_A,editTp1,(centreFreq-editOnFreq2).*Bfield.*gamma/1e6); %OFF
    editRFonB       =rf_freqshift(editRF_A,editTp1,(centreFreq-editOnFreq1).*Bfield.*gamma/1e6); %1.90ppm
    
    % Load the spin system definitions and pick the spin system of choice
    load spinSystems
    sys=eval(['sys' spinSys]);
    
    % Resample the refocusing RF pulses to 100 pts to reduce computational workload
    %   excRF = rf_resample(excRF,100);
      refRF = rf_resample(refRF,100);
    
    % Set up gradients
    
    % Gradient pulse shape generated from mpqspy (PPE) and divided by 10 to
    % convert mT/m to G/cm, scnh and mgs
    
    %Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    %Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
    
    %% Built from sim_MP_Siemens.m and sim_sLASER_shaped.m
    %Initialize structures:
    
    outA_posx_rpc = cell(length(x));
    outB_posx_rpc = cell(length(x));
    d_A           = struct([]);
    d_B           = struct([]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %loop through space: Don't forget to initialize the parallel processing
    %toolbox workers using 'matlabpool open N' (for N workers, 12 max).
    
    %if you do not have the parallel computing toolbox, uncomment the first
    %for loop and delete "parfor X=1:length(x)"
    %for X=1:length(x)
    parfor X=1:length(x)  %Use this if you do have the MATLAB parallel processing toolbox
%         for EP1=1:length(editPhCyc1)
%             for m=1:length(ph1)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x))]);
%                 ', ' ...
%                     'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ...
%                     '; Adiabatic phase cycle ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
                
                outA_posx_rpc{X} =sim_sLASER_shaped_fastX_fastphase(Bfield,sys,refRF,refTp,editRFonA,editTp1,x(X),Gx,editPhCyc1,ph1,ph2,taus,centreFreq,flipAngle);
                outB_posx_rpc{X} =sim_sLASER_shaped_fastX_fastphase(Bfield,sys,refRF,refTp,editRFonB,editTp1,x(X),Gx,editPhCyc1,ph1,ph2,taus,centreFreq,flipAngle);
                
%             end %end of 1st refocusing phase cycle loop.
%         end %end of 1st editing phase cycle loop.
    end %end of spatial loop (parfor) in x direction.
    
    %calculate the average density matrix (Doing this inside a separate for
    %loop because I couldn't figure out how to do this inside the parfor loop):
    for X=1:length(x)
%         for EP1=1:length(editPhCyc1)
%             for m=1:length(ph1)
               
                d_A=sim_dAdd(d_A,outA_posx_rpc{X});
                d_B=sim_dAdd(d_B,outB_posx_rpc{X});
            
%             end
%         end
    end
    
    outAA_posxy_rpc=cell(length(y));
    outBB_posxy_rpc=cell(length(y));
    outA           = struct([]);
    outB           = struct([]);
        
    %for Y=1:length(y)
    parfor Y=1:length(y)
        %         for EP1=1:length(editPhCyc1)
        %             for EP2=1:length(editPhCyc2)
        %                 for m=1:length(ph1)
        disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y))]);
        %                     ', '...
        %                         'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
        %                         'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
        %                         'Adiabatic phase cycle ' num2str(m) ' of ' num2str(length(ph1)) '!!!']);
        %                                                     sim_sLASER_shaped_fastYY(d,            n,sw,Bfield,linewidth,sys,te,RF,   refTp,editPulse,editTp,dy,   Gy,    editPh2,       ph3,ph4,taus,centreFreq,flipAngle)
        outAA_posxy_rpc{Y}=sim_sLASER_shaped_fastY_fastphase(d_A,Npts,sw,Bfield,lw,sys,TE,refRF,refTp,editRFonA,editTp1,y(Y),Gy,editPhCyc2,ph3,ph4,taus,centreFreq,flipAngle);
        outBB_posxy_rpc{Y}=sim_sLASER_shaped_fastY_fastphase(d_B,Npts,sw,Bfield,lw,sys,TE,refRF,refTp,editRFonB,editTp1,y(Y),Gy,editPhCyc2,ph3,ph4,taus,centreFreq,flipAngle);
        %                 end
        %             end
        %         end %end of 2nd editing phase cycle loop.
    end %end of 1st editing phase cycle loop.
    
    
    for Y=1:length(y)
        %         for EP1=1:length(editPhCyc1)
        %             for EP2=1:length(editPhCyc2)
        %                 for m=1:length(ph1)
        
        outA=op_addScans(outA,outAA_posxy_rpc{Y});
        outB=op_addScans(outB,outBB_posxy_rpc{Y});
        
        %                 end
        %             end
        %         end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For consistent scaling across different shaped simulations, we need to :
    %1.  Scale down by the total number of simulations run (since these were
    %    all added together.
    numSims=(nX*nY*length(editPhCyc1)*length(editPhCyc2)*length(ph1)*length(ph2)*length(ph3)*length(ph4));
    outA = op_ampScale(outA,1/numSims);
    outB = op_ampScale(outB,1/numSims);
    
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    if fovX>thkX
        voxRatio=(thkX*thkY)/(fovX*fovY);
    else
        voxRatio=1;
    end
    %out=op_ampScale(out,1/voxRatio);
    outA = op_ampScale(outA,1/voxRatio);
    outB = op_ampScale(outB,1/voxRatio);
    
    %Correct residual DC offset
    outA = op_dccorr(outA,'p');
    outB = op_dccorr(outB,'p');
    
    outA.name = metabolite;
    outB.name = metabolite;
    outA.nX=nX;
    outB.nX=nX;
    outA.thkX=thkX;
    outB.thkX=thkX;
    
    save(out_name,'outA','outB');
    
end
end %END OF MAIN FUNCTION:  Nested functions below.

%% In X direction
function d = sim_sLASER_shaped_fastX_fastphase(Bfield,sys,RF,refTp,editPulse,editTp,dx,Gx,editPh1,ph1,ph2,taus,centreFreq,flipAngle)

% if nargin<23
%flipAngle=180;
%     if nargin<22
%centreFreq=3.0; % Change it to GABA 3 ppm -- SH and MGSaleh
%     end
% end

%Check if this is a gradient modulated pulse.  If so, set Gx equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gx and then set Gx equal to zero:
    RF=rf_scaleGrad(RF,Gx);
    Gx=0;
end

%editFlip{1}=[0 0 180 180 0 0]; %Cell array of edit-on pulse flip angles.

delays=zeros(size(taus));
delays(1)=taus(1)-(editTp/2);
delays(2)=taus(2)-((editTp+refTp)/2);
delays(3)=taus(3)-((refTp+refTp)/2);
delays(4)=taus(4)-((refTp+refTp)/2);
delays(5)=taus(5)-((refTp+editTp)/2);
delays(6)=taus(6)-((editTp+refTp)/2);
delays(7)=taus(7)-(refTp/2);

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

% %BEGIN sLASER PULSE SEQUENCE************  Shaped editing pulse
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d_in=sim_evolve(d,H,delays(1)/1000);                       %Evolve by tau1
%% start of the editing pulse
clear d_out1;

for kk=1:length(editPh1)
    d_out1{kk}=sim_shapedRF(d_in,H,editPulse,editTp,flipAngle,editPh1(kk));          
end

for RP1=1:length(editPh1)
    if RP1==1
        out_edit1=d_out1{RP1};
    else
        out_edit1=sim_dAdd(out_edit1,d_out1{RP1});
    end
    
end
% end of 1st editing pulse

d_in=sim_evolve(out_edit1,H,delays(2)/1000);   %Evolve by tau2

%% 1st adiabatic refocusing pulse
for kk=1:length(ph1)
    d_out1{kk}=sim_shapedRF(d_in,H,RF,refTp,flipAngle,ph1(kk),dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
    
end

for RP1=1:length(ph1)
    receiver_phase=round((exp(1i*ph1(RP1)/180*pi*2)));
    %out_1=sim_dAdd(out_1,d_out1{RP1},receiver_phase);
    
    if RP1==1
        out_1=d_out1{RP1};
    else
        out_1=sim_dAdd(out_1,d_out1{RP1},receiver_phase);
    end

end
% end of 1st adiabatic refocusing pulse

d_in=sim_evolve(out_1,H,delays(3)/1000);                       %Evolve by tau2
%% 2nd adiabatic refocusing pulse
%clear d_in;
for kk=1:length(ph2)
    d_out2{kk}=sim_shapedRF(d_in,H,RF,refTp,flipAngle,ph2(kk),dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
end

for RP2=1:length(ph2)
    receiver_phase=round((exp(1i*ph2(RP2)/180*pi*2)));
    if RP2==1
        out_2=d_out2{RP2};
    else
        out_2=sim_dAdd(out_2,d_out2{RP2},receiver_phase);
    end    
end
% end of the editing refocusing pulse
d=sim_evolve(out_2,H,delays(4)/1000);                       %Evolve by tau2
%END PULSE SEQUENCE**************

%After running this many times along x, the density matrices should be
%averaged, and then the average density matrix should be passed through
%'sim_sLASER_shaped_fastY' at various different y-positions.



end

%% In Y direction
function out = sim_sLASER_shaped_fastY_fastphase(d,n,sw,Bfield,linewidth,sys,te,RF,refTp,editPulse,editTp,dy,Gy,editPh2,ph3,ph4,taus,centreFreq,flipAngle)

% if nargin<24
%flipAngle=180;
%     if nargin<23
%centreFreq=3.0; % Change it to GABA 3 ppm -- SH and MGSaleh
%     end
% end

if RF.isGM
    %Scale the GM waveform by the factor Gy and then set Gy equal to zero:
    RF=rf_scaleGrad(RF,Gy);
    Gy=0;
end

delays=zeros(size(taus));
delays(1)=taus(1)-(editTp/2);
delays(2)=taus(2)-((editTp+refTp)/2);
delays(3)=taus(3)-((refTp+refTp)/2);
delays(4)=taus(4)-((refTp+refTp)/2);
delays(5)=taus(5)-((refTp+editTp)/2);
delays(6)=taus(6)-((editTp+refTp)/2);
delays(7)=taus(7)-(refTp/2);

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end
%
%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

% %BEGIN sLASER PULSE SEQUENCE************  Shaped editing pulse
%% 3rd adiabatic refocusing pulse along Y gradient
for kk=1:length(ph3)
    d_out3{kk}=sim_shapedRF(d,H,RF,refTp,flipAngle,ph3(kk),dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
end
        for RP3=1:length(ph3)
            receiver_phase=round((exp(1i*ph3(RP3)/180*pi*2)));
             if RP3==1
                 out_3=d_out3{RP3};
             else
                  out_3=sim_dAdd(out_3,d_out3{RP3},receiver_phase);
             end
        end
        
d=sim_evolve(out_3,H,delays(5)/1000);                       %Evolve by tau2

%% 2nd editing pulse along Y gradient
for kk=1:length(editPh2) 
d_edit2{kk}=sim_shapedRF(d,H,editPulse,editTp,flipAngle,editPh2(kk));          
end

for RP4=1:length(editPh2)
    
    if RP4==1
        out_edit2=d_edit2{RP4};
    else
        out_edit2=sim_dAdd(out_edit2,d_edit2{RP4});
        
    end
end
d=sim_evolve(out_edit2,H,delays(6)/1000);                       %Evolve by tau1
%% 4th adiabatic refocusing pulse along Y gradient
for kk=1:length(ph4) 
d_out4{kk}=sim_shapedRF(d,H,RF,refTp,flipAngle,ph4(kk),dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
end

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

d=sim_evolve(out_4,H,delays(7)/1000);                       %Evolve by tau1
[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.68-centreFreq);

%Fill in structure header fields:
out.seq='MEGA_sLASER_FF';
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

