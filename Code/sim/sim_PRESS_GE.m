% simMP.m
% Georg Oeltzschner, Johns Hopkins University 2018
% adapted from run_simHERMES_GABAGSHCosmod_basis_sets_GABA_GLU_GLN_NAA
% Diana Rotaru, King's College London 2017
% adapted from
% run_simMegaPressShaped.m
% Jamie Near, McGill University 2014.
%
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment using GE product sequence with fully shaped editing
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
% editWaveform      = name of editing pulse waveform.
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

function [out] = sim_PRESS_GE(metabolite)

%addpath(genpath('/Users/Georg/Documents/MATLAB/FID-A'));

for kk=1:numel(metabolite)
    metabolite = metabolite{kk};
    
    % ********************PARAMETERS**********************************
    
    % Spin system to simulate
    spinSys     = metabolite;
    
    % General properties of the simulations
    Npts    = 8192;     % number of spectral points
    sw      = 4000;     % spectral width [Hz]
    Bfield  = 3;  % magnetic field strength [Tesla]
    lw      = 2;        % linewidth of the output spectrum [Hz]
    gamma   = 42577000; % gyromagnetic ratio
    TE      = 30; % TE [ms]
    tau1    = 15; %TE1 for first spin echo [ms]
    tau2    = TE - tau1; %TE2 for second spin echo [ms]
    
    % Output name
    out_name    = ['GE_PRESS_' num2str(TE) '_' spinSys '.mat'];
    
    % Define the pulse waveforms here
    exciteWaveform      = 'lpe90_ws_8520.pta';       % name of excitation pulse wavefo
    refocWaveform       = 'rfa_3.9ms.pta';         % name of refocusing pulse waveform.
    
    % Define frequency parameters for editing targets
    centerFreq  = 3;                 % Center frequency of MR spectrum [ppm] - water signal, stays the same on Siemens
    
    % Define pulse durations and flip angles
    excTp       = 3.6;              % duration of excitation pulse [ms], set to zero for hard pulse
    iso_delay   = 0.5123*excTp;    % time from phase centre to end of excitation pulse ("asymmetry factor") [ms]
    refTp       = 5.2;                % duration of refocusing pulses [ms], set to zero for hard pulse
    flips       = [137,137];        % flip angles of first an second refocusing pulses [degrees]
    
    % Define phase cycling pattern
    refPhCyc1   = [0,90];           % phase cycling steps for 1st refocusing pulse [degrees]
    refPhCyc2   = [0,90];           % phase cycling steps for 2nd refocusing pulse [degrees]
    
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
        y=linspace(-fovY/2,fovY/2,nY); %X positions to simulate [cm]
    else
        y=0;
    end
    if nZ>1
        z=linspace(-fovZ/2,fovZ/2,nZ); %X positions to simulate [cm]
    else
        z=0;
    end
    
    %Load RF waveform
    refRF=io_loadRFwaveform(refocWaveform,'ref',0);
    
    % Define the amplitudes of the slice-selective gradients
    % Read: x gradient for a voxel with 3 cm length in x direction is
    % 0.105058 G/cm as specified in the sequence
    Gx  = 0.106*(3/thkX);   % Gradient during 1st refocusing pulse [G/cm]
    Gy  = 0.106*(3/thkY);   % Gradient during 2nd refocusing pulse [G/cm]
    Gz  = 0.17826*(3/thkZ);   % Gradient during excitation pulse [G/cm]
    
%     Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
%     Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
%     
    %Load spin systems
    load spinSystems
    sys=eval(['sys' spinSys]);
    
    %Resample refocusing RF pulse from 400 pts to 100 pts to reduce
    %computational workload
    refRF=rf_resample(refRF,100);

    %Initialize structures:
    d_temp=cell(length(x),length(refPhCyc1));
    d=cell(length(refPhCyc1));
    
    %loop through space: If you are using the parfor loops below, and you are
    %using an older version of MATLAB (e.g.R2012), don't forget to initialize
    %the parallel processing toolbox workers using 'matlabpool open N' (for N
    %workers, 12 max).  I don't think this is necessary for newer version of
    %MATLAB.
    
    %First loop through all x-positions, simulating only the first refocusing
    %pulse.
    %First loop through x direction (first refoc pulse only);
    
    %for RP1=1:length(refPhCycl1)  %Use this if you don't have the MATLAB parallel processing toolbox
    parfor X=1:length(x)  %Use this if you have the MATLAB parallel processing toolbox
        for RP1=1:length(refPhCyc1)
            disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) '!!!']);
            d_temp{X}{RP1}=sim_press_shaped_fastRef1(Bfield,sys,tau1,tau2,refRF,refTp,x(X),Gx,refPhCyc1(RP1),centerFreq,flips(1));
        end
    end
    
    %calculate the average density matrix (Doing this inside a separate for
    %loop because I couldn't figure out how to do this inside the parfor loop):
    for X=1:length(x)
        for RP1=1:length(refPhCyc1)
            d{RP1}=sim_dAdd(d{RP1},d_temp{X}{RP1});
        end
    end
    
    
    % %Initialize structures:
    out_temp=cell(length(y),length(refPhCyc1),length(refPhCyc2));
    out=struct([]);
    
    %Now loop through y direction (second refoc pulse only);
    %for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
    parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
        for RP1=1:length(refPhCyc1)
            for RP2=1:length(refPhCyc2)
                disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                    'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                    'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
                out_temp{Y}{RP1}{RP2}=sim_press_shaped_fastRef2(d{RP1},Npts,sw,Bfield,lw,sys,tau1,tau2,...
                    refRF,refTp,y(Y),Gy,refPhCyc2(RP2),centerFreq,flips(2));
            end
        end
    end
    
    %Now combine the outputs;  Again, doing this inside a separate for loop
    %becuase I can't figure out how to do this inside the parfor loop:
    for Y=1:length(y)
        for RP1=1:length(refPhCyc1)
            for RP2=1:length(refPhCyc2)
                out=op_addScans(out,out_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
            end
        end
    end
    
    %For consistent scaling across different shaped simulations, we need to :
    %1.  Scale down by the total number of simulations run (since these were
    %    all added together.
    numSims=(nX*nY*length(refPhCyc1)*length(refPhCyc2));
    out=op_ampScale(out,1/numSims);
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    voxRatio=(thkX*thkY)/(fovX*fovY);
    out=op_ampScale(out,1/voxRatio);
    out = op_dccorr(out,'p');
    % Correct residual DC offset
    
    out.name=metabolite;
    out.centerFreq = centerFreq;          
    save(out_name,'out');
end
end







%Nested Function #1
    function d = sim_press_shaped_fastRef1(Bfield,sys,tau1,tau2,RF,tp,dx,Gx,phCyc1,centerFreq,flipAngle)
        %
        % USAGE:
        % d = sim_press_shaped_fastRef1(n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dx,Gx,phCyc1,flipAngle)
        %
        % DESCRIPTION:
        % This function simulates only the first bit of the PRESS experiment, up to
        % the beginning of the second refocusing pulse.  The excitation is
        % simulated as an instantaneous rotation, and the refocusing pulse is
        % simulated as a shaped rotation.
        %
        % This code is designed to be used in highly-accelerated shaped simulations,
        % using the method described by Yan Zhang et al. Med Phys 2017;44(8):
        % 4169-78.
        %
        % This code enables the choice of the phase of the refocusing pulse.  This
        % enables phase cycling of the refocusing pulses by repeating simulations
        % with different editing pulse phases, which is necessary to remove phase
        % artefacts from the editing pulses.  A four step phase cycling scheme is typically
        % sufficient, where both refocusing pulses are phase cycled by 0 and 90 degrees, and
        % the phase are combined in the following way:
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
        % n         = number of points in fid/spectrum
        % sw        = desired spectral width in [Hz]
        % Bfield    = main magnetic field strength in [T]
        % linewidth = linewidth in [Hz]
        % sys       = spin system definition structure
        % tau1      = echo time 1 in [ms].
        % tau2      = echo time 2 in [ms].
        % RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
        % tp        = RF pulse duration in [ms]
        % dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
        % dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
        % Gx        = gradient strength for first selective refocusing pulse [G/cm]
        % Gy        = gradient strength for second selective refocusing pulse [G/cm]
        % phCycl    = initial phase of the first refocusing pulse in [degrees];
        % phCycl2   = initial phase of the second refocusing pulse in [degrees];
        % flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
        %
        % OUTPUTS:
        % out       = simulated spectrum, in FID-A structure format, using PRESS
        %             sequence.
        
        if nargin<11
            flipAngle=180;
        end
        
        if tau1<tp/1000
            error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
        end
        if tau2<tp/1000
            error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
        end
        
        %Set water to centre
        centreFreq=centerFreq;
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
        
        %BEGIN PULSE SEQUENCE************
        d=sim_excite(d,H,'x');                                    %EXCITE
        d=sim_evolve(d,H,delays(1)/2000);                            %Evolve by delays(1)/2
        d=sim_shapedRF(d,H,RF,tp,flipAngle,90+phCyc1,dx,Gx);          %1st shaped 180 degree refocusing pulse
        d=sim_evolve(d,H,(delays(1)+delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
        %END PULSE SEQUENCE**************
        
        %After running this many times along x, the density matrices should be
        %averaged, and then the average density matrix should be passed through
        %'sim_press_shaped_fastRef2' at various different y-positions.
        
        
    end









%Nested Function #2
    function out = sim_press_shaped_fastRef2(d,n,sw,Bfield,linewidth,sys,tau1,tau2,RF,tp,dy,Gy,phCyc2,centerFreq,flipAngle)
        %
        % USAGE:
        % out = sim_press_shaped_fastRef2(d,n,sw,Bfield,linewidth,sys,tau2,RF,tp,dy,Gy,phCyc2,centerFreq,flipAngle)
        %
        % DESCRIPTION:
        % This function simulates only the last bit of the PRESS experiment, from the
        % the beginning of the second refocusing pulse, to the end.  The refocusing
        %pulse is simulated as a shaped rotation.
        %
        % This code is designed to be used in highly-accelerated shaped simulations,
        % using the method described by Yan Zhang et al. Med Phys 2017;44(8):
        % 4169-78.
        %
        % This code enables the choice of the phase of the refocusing pulse.  This
        % enables phase cycling of the refocusing pulses by repeating simulations
        % with different editing pulse phases, which is necessary to remove phase
        % artefacts from the editing pulses.  A four step phase cycling scheme is typically
        % sufficient, where both refocusing pulses are phase cycled by 0 and 90 degrees, and
        % the phase are combined in the following way:
        %
        % signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
        %
        % where, in [X Y], X is the phase of the first refocusing pulse and Y is
        % the phase of the second refocusing pulse
        %
        % Finally, this code simulates the spectrum at a given point in space (y),
        % given the values of the slice selection gradient (Gy).  In order
        % to fully simulate the MEGA-PRESS experiment, you have to first run
        % sim_press_shaped_fastRef1.m at all points in space (x), followed by
        % this code, at all points in space (y).
        %
        % INPUTS:
        % d         = starting density matrix (obtained using 'sim_press_shaped_fastRef1.m')
        % n         = number of points in fid/spectrum
        % sw        = desired spectral width in [Hz]
        % Bfield    = main magnetic field strength in [T]
        % linewidth = linewidth in [Hz]
        % sys       = spin system definition structure
        % tau1      = echo time 1 in [ms].
        % tau2      = echo time 2 in [ms].
        % RF        = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
        % tp        = RF pulse duration in [ms]
        % dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
        % dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
        % Gx        = gradient strength for first selective refocusing pulse [G/cm]
        % Gy        = gradient strength for second selective refocusing pulse [G/cm]
        % phCycl    = initial phase of the first refocusing pulse in [degrees];
        % phCycl2   = initial phase of the second refocusing pulse in [degrees];
        % flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
        %
        % OUTPUTS:
        % out       = simulated spectrum, in FID-A structure format, using PRESS
        %             sequence.
        
        if nargin<16
            flipAngle=180;
        end
        
        if tau1<tp/1000
            error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
        end
        if tau2<tp/1000
            error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
        end
        
        %Set water to centre
        centreFreq=centerFreq;
        for k=1:length(sys)
            sys(k).shifts=sys(k).shifts-centreFreq;
        end
        
        %Calculate Hamiltonian matrices and starting density matrix.
        [H]=sim_Hamiltonian(sys,Bfield);
        
        %Calculate new delays by subtracting the pulse duration from tau1 and tau2;
        delays=zeros(2);
        delays(1)=tau1-tp;
        delays(2)=tau2-tp;
        if sum(delays<0)
            error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
        end
        
        %BEGIN PULSE SEQUENCE************
        d=sim_shapedRF(d,H,RF,tp,flipAngle,90+phCyc2,dy,Gy);          %2nd shaped 180 degree refocusing pulse
        d=sim_evolve(d,H,delays(2)/2000);                            %Evolve by delays(2)/2
        [out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
        %END PULSE SEQUENCE**************
        
        %Correct the ppm scale:
        out.ppm=out.ppm-(4.68-centreFreq);
        
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
    end