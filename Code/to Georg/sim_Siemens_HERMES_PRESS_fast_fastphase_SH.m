% sim_HERMES_orig_Philips_SH.m

% Written by Steve CN Hui and Muhammad Saleh, Johns Hopkins University, 2019
% adapted from sim_MP_Siemens.m and sim_hercules_1_SH.m
% Georg Oeltzschner, Johns Hopkins University, 2018
% adapted from run_simSLaserShaped.m
% Dana Goerzen, Jamie Near, McGill University, 2019
%
% DESCRIPTION:
% This script simulates a HERMES experiment using Philips product sequence with fully shaped editing
% and refocusing pulses.  Phase cycling of both the editing and refocusing
% pulses is performed.  Furthermore, simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations,
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations. To enable the use
% of the MATLAB parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
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

function [outA,outB,outC,outD] = sim_Siemens_HERMES_PRESS_fast_fastphase_SH(metabolites)

addpath(genpath('/Users/steve/Documents/MATLAB/FID-A'));
% addpath(genpath('/home/shui5/matlab/FID-A')); %This is for KKI - SH 07252019

for kk=1:numel(metabolites)
    metabolite = metabolites{kk};
    
    % ********************PARAMETERS**********************************
    
    % Spin system to simulate
    spinSys     = metabolite;
    out_name    = ['Siemens_HERMES_PRESS_41pts_' spinSys '.mat']; 
    
    % General properties of the simulations
    Npts    = 8192;     % number of spectral points
    sw      = 4000;     % spectral width [Hz]
    Bfield  = 2.89;     % Philips (3T)/Siemens (2.89T) magnetic field strength [Tesla]
    lw      = 1;        % linewidth of the output spectrum [Hz]
    gamma   = 42577000; % gyromagnetic ratio
    
    % Define the pulse waveforms here
    exciteWaveform      = 'univ_spreddenrex.pta';       % name of excitation pulse waveform.
    refocWaveform       = 'univ_eddenrefo.pta';         % name of refocusing pulse waveform.
    editWaveform1       = 'sl_univ_pulse.pta';  % name of 1st single editing pulse waveform. [4.56ppm]
    editWaveform2       = 'sl_univ_pulse.pta';  % name of 2nd single editing pulse waveform. [1.9ppm]
    editWaveform3       = 'dl_Siemens_4_56_1_90.pta';   % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
    editWaveform4       = 'sl_univ_pulse.pta';  % name of non-editing pulse waveform. [non-editing]
    
    % Define frequency parameters for editing targets
    flipAngle = 180;
    centreFreq  = 3.0;              % Center frequency of MR spectrum [ppm]
    editOnFreq1 = 4.56;             % Center frequency of 1st HERMES experiment [ppm]
    editOnFreq2 = 1.9;              % Center frequency of 2nd HERMES experiment [ppm]
    editOnFreq3 = (4.56+1.9)/2;     % Center frequency of 3rd HERMES experiment [ppm]
    editOnFreq4 = 10.0;             % Center frequency of 4th HERMES experiment [ppm]
       
    % Define pulse durations and flip angles specific for every vendor
    excTp       = 7.1296;        % duration of excitation pulse [ms], set to zero for hard pulse
    iso_delay   = 0;        % time from phase centre to end of excitation pulse ("asymmetry factor") [ms]
    refTp       = 6.8944;           % duration of refocusing pulses [ms], set to zero for hard pulse
    editTp1     = 20;               % duration of 1st editing pulse [ms]
    editTp2     = 20;   	        % duration of 2nd editing pulse [ms]
    editTp3     = 20;               % duration of 3rd editing pulse [ms]
    editTp4     = 20;               % duration of 4th editing pulse [ms]
    TE          = 80;               % Echo time [ms]
    TE1         = 6.55*2;            % TE1 [ms] (Use 6.96*2 for Philips Original and 6.55*2 for Universial/Siemens)
    TE2         = TE - TE1;         % TE2 [ms]
    flips       = [180,180];        % flip angles of first an second refocusing pulses [degrees]
    
    % Define phase cycling pattern - stays the same, there are two editing
    % pulses in MEGA-PRESS with the same frequency, different time
%     editPhCyc1  = [0 90];           % phase cycling steps for 1st editing pulse [degrees]
%     editPhCyc2  = [0 90 180 270];   % phase cycling steps for 2nd editing pulse [degrees]
%     refPhCyc1   = [0,90];           % phase cycling steps for 1st refocusing pulse [degrees]
%     refPhCyc2   = [0,90];           % phase cycling steps for 2nd refocusing pulse [degrees]
    editPhCyc1  = [0 90 180 270];   % phase cycling steps for 1st editing pulse [degrees]
    editPhCyc2  = [0 90 180 270];   % phase cycling steps for 2nd editing pulse [degrees]
    refPhCyc1   = [0 90 180 270];   % phase cycling steps for 1st refocusing pulse [degrees]
    refPhCyc2   = [0 90 180 270];   % phase cycling steps for 2nd refocusing pulse [degrees]
    excPh       = 0;                % phase cycling steps for excitation pulse [degrees]
    
    % Define spatial resolution of simulation grid
    fovX    = 4.5; %size of the full simulation Field of View in the x-direction [cm]
    fovY    = 4.5; %size of the full simulation Field of View in the y-direction [cm]
    fovZ    = 4.5; %size of the full simulation Field of View in the y-direction [cm]
    thkX    = 3;    % slice thickness of x refocusing pulse [cm]
    thkY    = 3;    % slice thickness of y refocusing pulse [cm]
    thkZ    = 3;    % slice thickness of z excitation pulse [cm]
    nX      = 41;   % number of spatial points to simulate in x direction
    nY      = 41;   % number of spatial points to simulate in y direction
    nZ      = 41;   % number of spatial points to simulate in z direction
    if nX>1
        x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
    else
        x=0;
    end
    if nY>1
        y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]
    else
        y=0;
    end
    if nZ>1
        z=linspace(-fovZ/2,fovZ/2,nZ); %z positions to simulate [cm]
    else
        z=0;
    end
%     
%     % Define the amplitudes of the slice-selective gradients
%     % Read: x gradient for a voxel with 3 cm length in x direction is
%     % 0.105058 G/cm as specified in the sequence
%     Gx  = 0.106*(3/thkX);   % Gradient during 1st refocusing pulse [G/cm]
%     Gy  = 0.106*(3/thkY);   % Gradient during 2nd refocusing pulse [G/cm]
%     Gz  = 0.17826*(3/thkZ);   % Gradient during excitation pulse [G/cm] 
   

% set up exact timing - based on 'normal' pulse set on Philips 3T - SH 07252019
tausA = [TE1/2,...                           %middleEXC pulse to middleREFOC1
         ((TE/2 + TE1/2)/2)-TE1/2,...        %middleREFOC1 to middle 1st EDITING
         (TE/2+TE1/2)-((TE/2 + TE1/2)/2),... %middle 1st EDITING to middle REFOC2
         (3*TE/4 + TE1/4)-((TE/2+TE1/2)),... %middle REFOC2 to middle 2nd EDITING
         TE-(3*TE/4 + TE1/4)];               %middle 2nd EDITING to the start of readout  
     
%initialize evolution times (from sim_sLASER_shaped.m)
        %tau1=(te/4-tp)/2;
        %tau2=te/4-tp;
 
    % Below: Alternative for automatic calculation; double-check timings
    %     tausA=[TE1/2 - refTp/2,...
    %         TE2/4 - refTp/2- editTp/2,...
    %         TE2/4  + TE1/2 - refTp/2 - editTp/2,...
    %         TE2/4 - refTp/2 - editTp/2,...
    %         TE2/4 - editTp/2]; %timing of the pulse sequence of Exp A [ms]
    %     tausB=[TE1/2 - refTp/2,...
    %         TE2/4 - refTp/2- editTpB/2,...
    %         TE2/4  + TE1/2 - refTp/2 - editTpB/2,...
    %         TE2/4 - refTp/2 - editTpB/2,...
    %         TE2/4 - editTpB/2]; %timing of the pulse sequence of Exp B[ms]
    %     tausC=[TE1/2 - refTp/2,...
    %         TE2/4 - refTp/2- editTpC/2,...
    %         TE2/4  + TE1/2 - refTp/2 - editTpC/2,...
    %         TE2/4 - refTp/2 - editTpC/2,...
    %         TE2/4 - editTpC/2]; %timing of the pulse sequence of Exp C[ms]
    %     tausD=[TE1/2 - refTp/2,...
    %         TE2/4 - refTp/2- editTpD/2,...
    %         TE2/4  + TE1/2 - refTp/2 - editTpD/2,...
    %         TE2/4 - refTp/2 - editTpD/2,...
    %         TE2/4 - editTpD/2]; %timing of the pulse sequence of Exp D[ms]
    
    % ************END of PARAMETERS**********************************
    
    
    % ********************SET UP SIMULATION**********************************
    
    % Load RF waveforms for excitation and refocusing pulses
    
%     excRF = io_loadRFwaveform(exciteWaveform,'exc',0); %DJL added %Come back later - SH 07252019
    
    switch refocWaveform
%       case 'gtst1203_sp.pta'
%           refRF=io_loadRFwaveform_klc(refocWaveform,'ref',0);
%           refRF.tbw = 38.280;
%           refRF.tw1 = 9.932;
        case 'gtst1203_sp.pta'
            refRF=io_loadRFwaveform(refocWaveform,'ref',0);
            refRF.tbw = 1.357*6.8944; %BW99 (kHz) * dur (ms)
        case 'univ_eddenrefo.pta'
            refRF=io_loadRFwaveform(refocWaveform,'ref',0);
            refRF.tbw = 1.342*7; %BW99 (kHz) * dur (ms)
            refTp     = 7.0; 
    end    
     
    % Load RF waveforms for editing pulses
%     editRF1     =   io_loadRFwaveform(editWaveform1,'inv',20);
%     editRF2     =   io_loadRFwaveform(editWaveform2,'inv',20);
%     editRF3     =   io_loadRFwaveform(editWaveform3,'inv',20);
%     editRF4     =   io_loadRFwaveform(editWaveform4,'inv',20);
%     Double the RF power for dual band pulses
%     editRF3.tbw = editRF1.tbw;
%     editRF3.tw1 = editRF1.tw1*2;
%     editRF4.tbw = editRF1.tbw;
%     editRF4.tw1 = editRF1.tw1*2;
    editRF_A      = io_loadRFwaveform(editWaveform1,'inv',0); % sg100_100
%    editRF_A.tw1  = 0.93;
    editRF_B      = io_loadRFwaveform(editWaveform2,'inv',0); % sg100_100
%    editRF_B.tw1  = 0.93;
    editRF_C      = io_loadRFwaveform(editWaveform3,'inv',0); % dl_Philips_4_56_1_90
    editRF_C.tw1  = 1.8107;
    editRF_D      = io_loadRFwaveform(editWaveform4,'inv',0); % sg100_100
%    editRF_D.tw1  = 0.93;
        
    % Construct the editing pulses from the waveforms and defined
    % frequencies
    editRFonA     = rf_freqshift(editRF_A,editTp1,(centreFreq-editOnFreq1)*Bfield*gamma/1e6); %4.56
    editRFonB     = rf_freqshift(editRF_B,editTp2,(centreFreq-editOnFreq2)*Bfield*gamma/1e6); %1.9
    editRFonC     = rf_freqshift(editRF_C,editTp3,(centreFreq-editOnFreq3)*Bfield*gamma/1e6); %(4.56 + 1.90)/2
    editRFonD     = rf_freqshift(editRF_D,editTp4,(centreFreq-editOnFreq4)*Bfield*gamma/1e6); %non-editing
    % HERCULES has the same editing pulse duration and timing for all sub-experiments:
    editTp = editTp1;
    taus   = tausA;
  %%
    % Load the spin system definitions and pick the spin system of choice
    load spinSystems
    sys=eval(['sys' spinSys]);
    
    % Resample the refocusing RF pulses to 100 pts to reduce computational workload 
%   excRF = rf_resample(excRF,100);
    refRF = rf_resample(refRF,100);
    
    % Set up gradients
    Gx = (refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy = (refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
        
    %Initialize structures:
    d_A_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_B_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_C_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_D_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_A=cell(length(editPhCyc1),length(refPhCyc1));
    d_B=cell(length(editPhCyc1),length(refPhCyc1));
    d_C=cell(length(editPhCyc1),length(refPhCyc1));
    d_D=cell(length(editPhCyc1),length(refPhCyc1));        

%% Copied from simHERCULES.m (Check with Georg whether the excitation (z dimension) necessary)

%     tic;
%     [ d,H,delays,iso_delay ] = sim_megapress_shapedall_varflip_step0(centreFreq,sys,excTp,taus,refTp,editTp,Bfield,iso_delay);
%     disp(['PREP done']);
%     
%     % STEP 1
%     % Apply excitation pulse + evolution for d1
%     step=1;
%     sum_all=cell(1,length(sys));
%     for kk=1:length(sys)
%         sum_all{kk} = zeros(size(d{kk}));
%     end
%     d1_temp=cell(length(z),1);
%     
%     parfor Z=1:length(z)
%         disp(['Z-pos ' num2str(Z) '/' num2str(length(z))]);
%         d1_temp{Z} = sim_megapress_shapedall_varflip_start(step,d,H,delays,Npts,sw,Bfield,lw,taus,...
%             flips,sys,centreFreq,editRFonA,editTp,editPhCyc1(1),editPhCyc1(2),...
%             refRF,refTp,excRF,excTp,iso_delay,excPh,Gz,Gx,Gy,x,y,z(Z),refPhCyc1(1),refPhCyc1(2));
%         
%     end
%     %calculate the average density matrix (Doing this inside a separate for 
%     %loop because I couldn't figure out how to do this inside the parfor loop):
%     for Z=1:length(z)
%         sum_all=sim_dAdd(sum_all,d1_temp{Z},1);
%     end
%     % scale back down
%     for kk=1:length(sys)
%         d1_avg{kk}=sum_all{kk}./length(z);
%     end

%     % STEP 2
%     % Apply 1st refo and edit pulses + evolution for d2 and d3
%     step=2;
%     sum_A=cell(1,length(sys));
%     for kk=1:length(sys)
%         sum_A{kk} = zeros(size(d{kk}));
%     end
%     sum_B=sum_A;
%     sum_C=sum_A;
%     sum_D=sum_A;

   %% Built from sim_MP_Siemens.m 
   
   parfor X=1:length(x)  %Use this if you do have the MATLAB parallel processing toolbox
       %for EP1=1:length(editPhCyc1)
       %for RP1=1:length(refPhCyc1)
       disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x))]);
       %                 ', '...
       %                     'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
       %                     'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) '!!!']);
       d_A_temp{X}=sim_megapress_shaped_fastRef1_fastphase(Bfield,taus,sys,editRFonA,editTp,editPhCyc1,refRF,refTp,Gx,x(X),refPhCyc1,flipAngle,centreFreq);
       d_B_temp{X}=sim_megapress_shaped_fastRef1_fastphase(Bfield,taus,sys,editRFonB,editTp,editPhCyc1,refRF,refTp,Gx,x(X),refPhCyc1,flipAngle,centreFreq);
       d_C_temp{X}=sim_megapress_shaped_fastRef1_fastphase(Bfield,taus,sys,editRFonC,editTp,editPhCyc1,refRF,refTp,Gx,x(X),refPhCyc1,flipAngle,centreFreq);
       d_D_temp{X}=sim_megapress_shaped_fastRef1_fastphase(Bfield,taus,sys,editRFonD,editTp,editPhCyc1,refRF,refTp,Gx,x(X),refPhCyc1,flipAngle,centreFreq);
       %end %end of 1st refocusing phase cycle loop.
       %end %end of 1st editing phase cycle loop.
   end %end of spatial loop (parfor) in x direction.
   
    %calculate the average density matrix (Doing this inside a separate for
    %loop because I couldn't figure out how to do this inside the parfor loop):
    for X=1:length(x)
        %for EP1=1:length(editPhCyc1)
        %for RP1=1:length(refPhCyc1)
        d_A=sim_dAdd(d_A,d_A_temp{X});
        d_B=sim_dAdd(d_B,d_B_temp{X});
        d_C=sim_dAdd(d_C,d_C_temp{X});
        d_D=sim_dAdd(d_D,d_D_temp{X});
        %end
        %end
    end
    
    
    % %Initialize structures:
    outA_temp=cell(length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
    outB_temp=cell(length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
    outC_temp=cell(length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
    outD_temp=cell(length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
    outA=struct([]);
    outB=struct([]);
    outC=struct([]);
    outD=struct([]);
    
    %Now loop through y direction (second refoc pulse only);
    %for Y=1:length(y)  %Use this if you don't have the MATLAB parallel processing toolbox
    parfor Y=1:length(y)  %Use this if you do have the MATLAB parallel processing toolbox
        %         for EP1=1:length(editPhCyc1)
        %             for EP2=1:length(editPhCyc2)
        %                 for RP1=1:length(refPhCyc1)
        %                     for RP2=1:length(refPhCyc2)
        disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y))]);
        %                         ', '...
        %                             'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
        %                             'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
        %                             'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
        %                             'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
        outA_temp{Y}=sim_megapress_shaped_fastRef2_fastphase(d_A,Npts,sw,Bfield,lw,taus,sys,editRFonA,editTp,editPhCyc2,refRF,refTp,Gy,y(Y),refPhCyc2,flipAngle,centreFreq);
        outB_temp{Y}=sim_megapress_shaped_fastRef2_fastphase(d_B,Npts,sw,Bfield,lw,taus,sys,editRFonB,editTp,editPhCyc2,refRF,refTp,Gy,y(Y),refPhCyc2,flipAngle,centreFreq);
        outC_temp{Y}=sim_megapress_shaped_fastRef2_fastphase(d_C,Npts,sw,Bfield,lw,taus,sys,editRFonC,editTp,editPhCyc2,refRF,refTp,Gy,y(Y),refPhCyc2,flipAngle,centreFreq);
        outD_temp{Y}=sim_megapress_shaped_fastRef2_fastphase(d_D,Npts,sw,Bfield,lw,taus,sys,editRFonD,editTp,editPhCyc2,refRF,refTp,Gy,y(Y),refPhCyc2,flipAngle,centreFreq);
        %                     end %end of 2nd refocusing phase cycle loop.
        %                 end %end of 1st refocusing phase cycle loop.
        %             end %end of 2nd editing phase cycle loop.
        %         end %end of 1st editing phase cycle loop.
    end %end of spatial loop (parfor) in y direction.
    
    %Now combine the outputs;  Again, doing this inside a separate for loop
    %becuase I can't figure out how to do this inside the parfor loop:
    for Y=1:length(y)
%         for EP1=1:length(editPhCyc1)
%             for EP2=1:length(editPhCyc2)
%                 for RP1=1:length(refPhCyc1)
%                     for RP2=1:length(refPhCyc2)
                        outA=op_addScans(outA,outA_temp{Y});
                        outB=op_addScans(outB,outB_temp{Y});
                        outC=op_addScans(outC,outC_temp{Y});
                        outD=op_addScans(outD,outD_temp{Y});
%                     end
%                 end
%             end
%         end
    end
    
    %For consistent scaling across different shaped simulations, we need to :
    %1.  Scale down by the total number of simulations run (since these were
    %    all added together.
    numSims=(nX*nY*length(refPhCyc1)*length(refPhCyc2)*length(editPhCyc1)*length(editPhCyc2));
    outA=op_ampScale(outA,1/numSims);
    outB=op_ampScale(outB,1/numSims);
    outC=op_ampScale(outC,1/numSims);
    outD=op_ampScale(outD,1/numSims);
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    if fovX>thkX
        voxRatio=(thkX*thkY)/(fovX*fovY);
    else
        voxRatio=1;
    end
    outA=op_ampScale(outA,1/voxRatio);
    outB=op_ampScale(outB,1/voxRatio);
    outC=op_ampScale(outC,1/voxRatio);
    outD=op_ampScale(outD,1/voxRatio);
    
    % Correct residual DC offset
    outA = op_dccorr(outA,'p');
    outB = op_dccorr(outB,'p');
    outC = op_dccorr(outC,'p');
    outD = op_dccorr(outD,'p');
    
    outA.name=metabolite;
    outB.name=metabolite;
    outC.name=metabolite;
    outD.name=metabolite;
    
    outA.centreFreq=centreFreq;
    outB.centreFreq=centreFreq;
    outC.centreFreq=centreFreq;
    outD.centreFreq=centreFreq;
     
    outA.nX=nX;
    outB.nX=nX;
    outC.nX=nX;
    outD.nX=nX;
    
    outA.thkX=thkX;
    outB.thkX=thkX;
    outC.thkX=thkX;
    outD.thkX=thkX;
    
    save(out_name,'outA','outB','outC','outD');
    
end
end %END OF MAIN FUNCTION:  Nested functions below.

%Nested Function #1
function d = sim_megapress_shaped_fastRef1_fastphase(Bfield,taus,sys,editPulse,editTp,editPh1,refPulse,refTp,Gx,dx,refPh1,flipAngle,centreFreq)
%
% USAGE:
% out = sim_megapress_shaped_fastRef1(Bfield,taus,sys,editPulse,editTp,editPh1,refPulse,refTp,Gx,dx,refPh1,flipAngle)
%
% DESCRIPTION:
% This function simulates only the first bit of the MEGA-PRESS experiment,
% up to the beginning of the second refocusing pulse.  The excitation is
% simulated as an instantaneous rotation, and the refocusing and editing
% pulses are simulated as shaped rotations.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8):
% 4169-78.
%
% This code enables the choice of the phase of the refocusing and editing
% pulses.  This enables phase cycling of the refocusing pulses by repeating
% simulations with different refocusing and editing pulse phases, which is
% necessary to remove phase artefacts from the editing pulses.  A four step
% phase cycling scheme is typically sufficient, where both refocusing
% pulses are phase cycled by 0 and 90 degrees, and the phase are combined
% in the following way:
%
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
%
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
%
% Finally, this code simulates the spectrum at a given point in space (x),
% given the values of the slice selection gradient (Gx).  In order
% to fully simulate the MEGA-PRESS experiment, you have to run this
% simulation many times at various points in space (x), followed by
% sim_press_shaped_fastRef2.m, at all points in space (y).
%
% INPUTS:
% Bfield    = main magnetic field strength in [T]
% taus(1)     = time in [ms] from 90 to 1st 180
% taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
% taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
% taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
% taus(5)     = time in [ms] from 2nd edit pulse to ADC
%               FOR MEGA-PRESS on SIEMENS SYSTEM:  taus=[4.545,12.7025,21.7975,12.7025,17.2526];
% sys        = Metabolite spin system definition structure;
% editPulse  = RF pulse definition structure for editing pulses (obtain using 'io_loadRFwaveform.m')
% editTp     = duration of editing pulse in [ms];
% editPh1    = the phase of the first editing pulse in [degrees];
% refPulse   = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% refTp      = duration of refocusing pulse in [ms]
% Gx         = gradient strength for first selective refocusing pulse [G/cm]
% dx         = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% refPh1     = the phase of the first refocusing pulse in [degrees];
% flipAngle  = the flip angle of the refocusing pulse.
%
% OUTPUTS:
% d          = output density matrix

% if nargin<12
%     flipAngle=180;
% end
% 
% %Set 3ppm GABA resonance to centre
% centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-((refTp+editTp)/2);
delays(3)=taus(3)-((editTp+refTp)/2);
delays(4)=taus(4)-((refTp+editTp)/2);
delays(5)=taus(5)-(editTp/2);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end


%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                          %EXCITE
d=sim_evolve(d,H,delays(1)/1000);                               %Evolve by delays(1)
%d=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh1,dx,Gx);   %1st shaped 180 degree refocusing pulse
clear d_out1;
%d_in=d;
for kk=1:length(refPh1)
    d_out{kk}=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh1(kk),dx,Gx);   %1st shaped 180 degree refocusing pulse
    % average as we go
end
receiver_phase = exp(1i*refPh1/180*pi*2);
%Copy from FID-A , calculate the average density matrix (Doing this inside a separate for loop): 
for RP1=1:length(refPh1)
   if RP1==1
       d_out1=d_out{RP1};
   else
       d_out1=sim_dAdd(d_out1,d_out{RP1},receiver_phase(RP1));
   end
end

d=sim_evolve(d_out1,H,delays(2)/1000);                               %Evolve by delays(2)
%d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh1);            %1st shaped editing pulse rotation
clear d_out1;

for kk=1:length(editPh1)
        d_out1{kk}=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh1(kk));            %1st shaped editing pulse rotation
end
% end loop over editph1 and average over editph1

%out_edit1=struct([]);

%Now combine the outputs;
for RP1=1:length(editPh1)
    %receiver_phase=round((exp(1i*editPh1(RP1)/180*pi*2)));
    %out_edit1=sim_dAdd(out_edit1,d_out1{RP1});
    if RP1==1
        out_edit1=d_out1{RP1};
    else
        out_edit1=sim_dAdd(out_edit1,d_out1{RP1});
    end
    
end

d=sim_evolve(out_edit1,H,delays(3)/1000);                               %Evolve by delays(3)
%END PULSE SEQUENCE**************


%After running this many times along x, the density matrices should be
%averaged, and then the average density matrix should be passed through
%'sim_megapress_shaped_fastRef2' at various different y-positions.


end



%Nested Function #2
function out = sim_megapress_shaped_fastRef2_fastphase(d,n,sw,Bfield,linewidth,taus,sys,editPulse,editTp,editPh2,refPulse,refTp,Gy,dy,refPh2,flipAngle,centreFreq)
%
% USAGE:
% out = sim_press_shaped_fastRef2(d,n,sw,Bfield,linewidth,taus,sys,editPulse,editTp,editPh2,refPulse,refTp,Gy,dy,refPh2,flipAngle)
%
% DESCRIPTION:
% This function takes a starting density matrix, and simulates only the
% last bit of the MEGA-PRESS experiment, from the second refocusing pulse
% to the end.  The refocusing and editing pulses are simulated as shaped
% rotations.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8):
% 4169-78.
%
% This code enables the choice of the phase of the refocusing and editing
% pulses.  This enables phase cycling of the refocusing pulses by repeating
% simulations with different refocusing and editing pulse phases, which is
% necessary to remove phase artefacts from the editing pulses.  A four step
% phase cycling scheme is typically sufficient, where both refocusing
% pulses are phase cycled by 0 and 90 degrees, and the phase are combined
% in the following way:
%
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
%
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
%
% Before running this code, you need to run sim_megapress_shaped_fast1.m to
% generate the density matrix resulting from the first part of the
% megapress sequence over a range of x positions.  Then you need to take
% the average density matrix across x.  Finally, this code simulates the
% spectrum at a given point in space (y), given the values of the slice
% selection gradient (Gy).  In order to fully simulate the MEGA-PRESS
% experiment, you have to run this simulation many times at various points
% in space (y).
%
% INPUTS:
% Bfield    = main magnetic field strength in [T]
% taus(1)     = time in [ms] from 90 to 1st 180
% taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
% taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
% taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
% taus(5)     = time in [ms] from 2nd edit pulse to ADC
%               FOR MEGA-PRESS on SIEMENS SYSTEM:  taus=[4.545,12.7025,21.7975,12.7025,17.2526];
% sys        = Metabolite spin system definition structure;
% editPulse  = RF pulse definition structure for editing pulses (obtain using 'io_loadRFwaveform.m')
% editTp     = duration of editing pulse in [ms];
% editPh1    = the phase of the first editing pulse in [degrees];
% refPulse   = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% refTp      = duration of refocusing pulse in [ms]
% Gx         = gradient strength for first selective refocusing pulse [G/cm]
% dx         = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% refPh1     = the phase of the first refocusing pulse in [degrees];
% flipAngle  = the flip angle of the refocusing pulse.
%
% OUTPUTS:
% d          = output density matrix

% if nargin<16
%     flipAngle=180;
% end
% 
% %Set 3ppm GABA resonance to centre
% centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-((refTp+editTp)/2);
delays(3)=taus(3)-((editTp+refTp)/2);
delays(4)=taus(4)-((refTp+editTp)/2);
delays(5)=taus(5)-(editTp/2);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end


%BEGIN PULSE SEQUENCE************
%d=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh2,dy,Gy);  %2nd shaped 180 degree refocusing pulse
clear d_in d_temp d_out1;
d_in=d;
for kk=1:length(refPh2)
    d_temp{kk}=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh2(kk),dy,Gy);  %2nd shaped 180 degree refocusing pulse
% average as we go
end
receiver_phase = exp(1i*refPh2/180*pi*2);
%Copy from FID-A , calculate the average density matrix (Doing this inside a separate for loop): 
for RP2=1:length(refPh2)
   if RP2==1
       d_out1=d_temp{RP2};
   else
       d_out1=sim_dAdd(d_out1,d_temp{RP2},receiver_phase(RP2));
   end
end

d=sim_evolve(d_out1,H,delays(4)/1000);                          %Evolve by delays(4)
%d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh2);     %2nd shaped editing pulse rotation
for kk=1:length(editPh2) 
d_edit2{kk}=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh2(kk));     %2nd shaped editing pulse rotation       
end

%out_edit2=struct([]);

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:

for RP2=1:length(editPh2)
    %receiver_phase=round((exp(1i*editPh2(RP4)/180*pi*2)));
    %out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
    %out_edit2=sim_dAdd(out_edit2,d_edit2{RP4});
    
    if RP2==1
        out_edit2=d_edit2{RP2};
    else
        out_edit2=sim_dAdd(out_edit2,d_edit2{RP2});
        
    end
end

d=sim_evolve(out_edit2,H,delays(5)/1000);                          %Evolve by delays(5)
[out,~]=sim_readout(d,H,n,sw,linewidth,90);           %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.68-centreFreq);

%Fill in structure header fields:
out.seq='HERMES_PRESS_FF';
out.te=sum(taus);
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

