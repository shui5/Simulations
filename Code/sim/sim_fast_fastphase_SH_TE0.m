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


function [outA] = sim_fast_fastphase_SH_TE0(metabolites)

addpath(genpath('/Users/steve/Documents/MATLAB/FID-A'));
% addpath(genpath('/home/shui5/matlab/FID-A')); %This is for KKI - SH 07252019

for kk=1:numel(metabolites)
    metabolite = metabolites{kk};
    
    % ********************PARAMETERS**********************************
    
    % Spin system to simulate
    spinSys     = metabolite;
    out_name    = ['TE_0_101pts_' spinSys '.mat']; 
    
    % General properties of the simulations
    Npts    = 8192;     % number of spectral points
    sw      = 4000;     % spectral width [Hz]
    Bfield  = 3;        % Philips magnetic field strength [Tesla]
    lw      = 1;        % linewidth of the output spectrum [Hz]
    gamma   = 42577000; % gyromagnetic ratio
    
    % Define the pulse waveforms here
    exciteWaveform      = 'univ_spreddenrex.pta';       % name of excitation pulse waveform.
    refocWaveform       = 'gtst1203_sp.pta';            % name of refocusing pulse waveform.
    editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [4.56ppm]
    editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 2nd single editing pulse waveform. [1.9ppm]
    editWaveform3       = 'dl_Philips_4_56_1_90.pta';   % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
    editWaveform4       = 'sg100_100_0_14ms_88hz.pta';  % name of non-editing pulse waveform. [non-editing]
    
    % Define frequency parameters for editing targets
    flipAngle = 180;
    centreFreq  = 3.0;              % Center frequency of MR spectrum [ppm]
    editOnFreq1 = 4.1;             % Center frequency of 1st HERMES experiment [ppm]
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
    TE          = 0;               % Echo time [ms]
    TE1         = 6.96*2            % TE1 [ms] (Use 6.96*2 for Philips Original and 6.55*2 for Universial/Siemens)
    TE2         = TE - TE1;         % TE2 [ms]
    flips       = [180,180];        % flip angles of first an second refocusing pulses [degrees]
    
    % Define phase cycling pattern - stays the same, there are two editing
    % pulses in MEGA-PRESS with the same frequency, different time

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
    nX      = 101;   % number of spatial points to simulate in x direction
    nY      = 101;   % number of spatial points to simulate in y direction
    nZ      = 101;   % number of spatial points to simulate in z direction
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
   

% set up exact timing - based on 'normal' pulse set on Philips 3T - SH 07252019
tausA = [TE1/2,...                           %middleEXC pulse to middleREFOC1
         ((TE/2 + TE1/2)/2)-TE1/2,...        %middleREFOC1 to middle 1st EDITING
         (TE/2+TE1/2)-((TE/2 + TE1/2)/2),... %middle 1st EDITING to middle REFOC2
         (3*TE/4 + TE1/4)-((TE/2+TE1/2)),... %middle REFOC2 to middle 2nd EDITING
         TE-(3*TE/4 + TE1/4)];               %middle 2nd EDITING to the start of readout  
    
    % ************END of PARAMETERS**********************************
    
    
    % ********************SET UP SIMULATION**********************************
    

  %%
    % Load the spin system definitions and pick the spin system of choice
    load spinSystems
    sys=eval(['sys' spinSys]);
            
    %Initialize structures:
    d_A_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_B_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_C_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_D_temp=cell(length(x),length(editPhCyc1),length(refPhCyc1));
    d_A=cell(length(editPhCyc1),length(refPhCyc1));
    d_B=cell(length(editPhCyc1),length(refPhCyc1));
    d_C=cell(length(editPhCyc1),length(refPhCyc1));
    d_D=cell(length(editPhCyc1),length(refPhCyc1));        

   %% Built from sim_MP_Siemens.m 
        outA=struct([]);

       outA=sim_megapress_shaped_fastRef1_fastphase(Bfield,sys,Npts,sw,lw,flipAngle,centreFreq);
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    if fovX>thkX
        voxRatio=1;
    else
        voxRatio=1;
    end

    
    % Correct residual DC offset
    outA = op_dccorr(outA,'p');

    outA.name=metabolite;
   
    outA.centreFreq=centreFreq;
    
    save(out_name,'outA'); % ,'outB','outC','outD');
    
end

x_lim = [1,5];

figure(1), plot(outA.ppm, outA.specs)
legend('Lac'),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('Lac TE0')

end %END OF MAIN FUNCTION:  Nested functions below.

%Nested Function #1
function out = sim_megapress_shaped_fastRef1_fastphase(Bfield,sys,n,sw,linewidth,flipAngle,centreFreq)

% %Set 3ppm GABA resonance to centre
% centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);


%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                          %EXCITE
[out,~]=sim_readout(d,H,n,sw,linewidth,90);           %Readout along y (90 degree phase); SCNH
%Correct the ppm scale:
out.ppm=out.ppm-(4.68-centreFreq);

%Fill in structure header fields:
out.seq='TE0';
out.te=0;
out.sim='ideal';

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
