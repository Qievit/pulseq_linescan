%% Line-scan DWI with TSE readout
% 
%
%
%

clear all
addpath('.\matlab\')

%% Initialize the parameters
outfile = '.\seqs\';
system = mr.opts('maxGrad', 45,'GradUnit', 'mT/m', ... % SKYRA
    'MaxSlew', 80, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 100e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 100e-6); 
seq=mr.Sequence(system);
gamma = 42.58e6/(2*pi) * 3; % [Hz] at 3T gyromagnetic ratio

feDir = 'x'; peDir = 'y'; ssDir = 'z';

% Sequence parameters
TE = 15e-3;
TR = 1;

fov = 192e-3; Nfe = 64; Nlin=1;
sliceThickness = 3e-3;
alphaRo = 180;
adcDur = 4.48e-3;
rfDur = 3e-3;

deltak = 1/fov;

%% DWI
dwType = 'mono'; % 'mono' 'dbipolar'
dwAmp = gamma*[ 40e-3 0 0;
                30e-3 0 0;
                40e-3 0 0]; % DWI gradient amplitude. FE PE SS [Hz/m] 
nb0 = 2;
dwDur = 12e-3; % DWI gradient duration. For bipolar: single lobe
Delta = 10e-3; % distance between diffusion gradients NOT including the gradient itself

dwAmp = [repmat([0 0 0],nb0,1) ; dwAmp];

%% Initialize RF pulses
[rfEx, gssExNorm] = mr.makeSincPulse(90*pi/180,system,'Duration',rfDur,'SliceThickness',sliceThickness, ...
    'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','excitation');

[rfRo, gssRo] = mr.makeSincPulse(alphaRo*pi/180,system,'Duration',rfDur,'SliceThickness',sliceThickness,...
    'apodization',0.5,'timeBwProduct',4,'use','refocusing');

rfRo.delay = system.rfDeadTime;

%% Initialize gradients, ADC, delays
% SS gradients
ssSpoilingArea = (4*pi)/(2*pi*sliceThickness); % 4 pi, gamma divides itself away
csMax = sqrt(ssSpoilingArea*system.maxSlew + gssRo.amplitude^2/2);
cs1 = ceil( csMax/system.maxSlew/system.gradRasterTime )*system.gradRasterTime;
cs2 = ceil( abs(csMax - gssRo.amplitude)/system.maxSlew/system.gradRasterTime )*system.gradRasterTime;
ssC1 = mr.makeExtendedTrapezoid(ssDir,system,'times',[0 cs1 cs1+cs2],'amplitudes',[0 csMax gssRo.amplitude]);
ssC2 = mr.makeExtendedTrapezoid(ssDir,system,'times',[0 cs2 cs1+cs2],'amplitudes',[gssRo.amplitude csMax 0]);
gssAdj = mr.makeExtendedTrapezoid(ssDir,system,'times',[0 rfDur+system.rfDeadTime+system.rfRingdownTime],'amplitudes',[gssRo.amplitude gssRo.amplitude]);

% FE gradients
gfeAmp = mr.makeTrapezoid(feDir,system,'FlatArea',Nfe*deltak,'FlatTime',adcDur);
feSpoilingArea = gfeAmp.area;
cfMaxA = system.maxGrad*0.5;
cf1 = ceil( cfMaxA/system.maxSlew/system.gradRasterTime )*system.gradRasterTime;
cf3 = ceil( (cfMaxA-gfeAmp.amplitude)/system.maxSlew/system.gradRasterTime )*system.gradRasterTime;
cf2 = ceil( ((feSpoilingArea - cf1*cfMaxA/2 - cf3*(cfMaxA-gfeAmp.amplitude)/2 - cf3*gfeAmp.amplitude)/cfMaxA)/system.gradRasterTime )*system.gradRasterTime;
feC1 = mr.makeExtendedTrapezoid(feDir,'times',[0 cf1 cf1+cf2 cf1+cf2+cf3],'amplitudes',[0 cfMaxA cfMaxA gfeAmp.amplitude]);
feC2 = mr.makeExtendedTrapezoid(feDir,'times',[0 cf3 cf2+cf3 cf1+cf2+cf3],'amplitudes',[gfeAmp.amplitude cfMaxA cfMaxA 0]);
gfe = mr.makeExtendedTrapezoid(feDir,'times',[0 adcDur+2*system.adcDeadTime],'amplitudes',[gfeAmp.amplitude gfeAmp.amplitude]);
adc = mr.makeAdc(Nfe,system,'Duration',adcDur,'Delay',system.adcDeadTime);

ieTE = rfDur + mr.calcDuration(gfe) + max([mr.calcDuration(feC1),mr.calcDuration(ssC1)])*2;
disp(['Using interecho spacing = ' num2str(ieTE*1e3) ' ms'])
delayRo = max([mr.calcDuration(feC1),mr.calcDuration(ssC1)]);

% Rephasing and prephasing
gfePreph = mr.makeTrapezoid(feDir,system,'Area',gfe.area/2 + feC1.area,'Duration',delayRo + mr.calcDuration(adc)/2);
gssReph = mr.makeTrapezoid(peDir,system,'Area',-gssExNorm.area/2,'Duration',delayRo + mr.calcDuration(adc)/2); % for the excitation pulse
gssPreph = mr.makeTrapezoid(ssDir,system,'Area',-(gssRo.area + ssSpoilingArea)/2,'Duration',delayRo + mr.calcDuration(adc)/2); % For the refocussing pulse
gssEx = mr.makeExtendedTrapezoid(peDir,'times',[0 gssExNorm.riseTime gssExNorm.riseTime+gssExNorm.flatTime mr.calcDuration(gssExNorm) mr.calcDuration(gssExNorm)+gssReph.riseTime mr.calcDuration(gssExNorm)+gssReph.riseTime+gssReph.flatTime mr.calcDuration(gssExNorm)+mr.calcDuration(gssReph)], ...
        'amplitudes',[0 gssExNorm.amplitude gssExNorm.amplitude 0 gssReph.amplitude gssReph.amplitude 0]);
gssPreph.delay = mr.calcDuration(gssExNorm);
gfePreph.delay = mr.calcDuration(gssExNorm);

% DW/ARFI gradients
if sum(strcmp(dwType,{'mono','sbipolar','dbipolar'}))
    for dir=1:size(dwAmp,1)
        dwTemplate = mr.makeTrapezoid(ssDir,system,'Amplitude',max(dwAmp(dir,:)),'Duration',dwDur/(1+sum(strcmp(dwType,{'sbipolar','dbipolar'}))));
        dwTimes = [0 dwTemplate.riseTime dwTemplate.riseTime+dwTemplate.flatTime mr.calcDuration(dwTemplate)];
        dwAmplitudes = [[0 0 0]' dwAmp(dir,:)' dwAmp(dir,:)' [0 0 0]'];
        if sum(strcmp(dwType,{'sbipolar','dbipolar'}))
            dwTimes = [dwTimes dwTimes(2:end)+mr.calcDuration(dwTemplate)];
            dwAmplitudes = [dwAmplitudes -dwAmplitudes(:,2:end)];
        end
        dwFe{dir} = mr.makeExtendedTrapezoid(feDir,system,'times',dwTimes,'amplitudes',dwAmplitudes(1,:));
        dwPe{dir} = mr.makeExtendedTrapezoid(peDir,system,'times',dwTimes,'amplitudes',dwAmplitudes(2,:));
        dwSs{dir} = mr.makeExtendedTrapezoid(ssDir,system,'times',dwTimes,'amplitudes',dwAmplitudes(3,:));
    end
end

feC1.delay = delayRo - mr.calcDuration(feC1);
TRdelay = TR ... 
         -(mr.calcDuration(gssEx,gssPreph,gfePreph) + Delta + dwFe{1}.shape_dur*2) ...
         -(delayRo*2 + mr.calcDuration(adc)) ...
         -(mr.calcDuration(gssAdj) + delayRo*2 + mr.calcDuration(adc))*Nlin;

%% Build sequence
seq.addBlock(mr.makeLabel('SET','SLC',0))
seq.addBlock(mr.makeLabel('SET','REP',0))

for dir = 1:size(dwAmp,1)
    seq.addBlock(rfEx,gssEx,gssPreph,gfePreph)
    
    switch dwType
        case 'mono'
            seq.addBlock(dwFe{dir},dwPe{dir},dwSs{dir})
            seq.addBlock(mr.makeDelay(Delta/2-mr.calcDuration(rfRo,gssRo)/2))
            seq.addBlock(rfRo,gssRo)
            seq.addBlock(mr.makeDelay(Delta/2-mr.calcDuration(rfRo,gssRo)/2))
            seq.addBlock(dwFe{dir},dwPe{dir},dwSs{dir})
            % seq.addBlock(mr.makeDelay(mr.calcDuration(ssC1)))
        case 'dbipolar'
            seq.addBlock(dwFe{dir},dwPe{dir},dwSs{dir})
            seq.addBlock(mr.makeDelay(Delta/2-mr.calcDuration(rfRo,gssAdj)/2-mr.calcDuration(ssC1)))
            seq.addBlock(ssC1)
            seq.addBlock(rfRo,gssAdj)
            seq.addBlock(ssC2)
            seq.addBlock(mr.makeDelay(Delta/2-mr.calcDuration(rfRo,gssAdj)/2-mr.calcDuration(ssC2)))
            seq.addBlock(dwFe{dir},dwPe{dir},dwSs{dir})
    end
    
        % ssC1.delay = mr.calcDuration(gfePreph) - mr.calcDuration(ssC1);
    % seq.addBlock(ssC1)
    %% Readout
    seq.addBlock(feC1,mr.makeLabel('SET', 'LIN', 1) )
    seq.addBlock(gfe,adc)
    ssC1.delay=0; ssC1.delay = delayRo - mr.calcDuration(ssC1);
    if Nlin-1 == 0; seq.addBlock(feC2);
    else;seq.addBlock(ssC1,feC2);end


    for i = 1:Nlin-1
        seq.addBlock(rfRo,gssAdj)
        seq.addBlock(ssC2,feC1,delayRo,mr.makeLabel('SET', 'LIN', 1) )
        seq.addBlock(gfe,adc)
        ssC1.delay=0; ssC1.delay = delayRo - mr.calcDuration(ssC1);
        if i==Nlin; seq.addBlock(feC2,mr.makeDelay(delayRo)); 
        else;seq.addBlock(ssC1,feC2,delayRo);end
    end
    min_seq_duration = sum(seq.blockDurations);
    seq.addBlock(mr.makeDelay(TRdelay))
    % seq.addBlock(mr.makeDelay(ceil((TR-min_seq_duration)/seq.gradRasterTime)*seq.gradRasterTime))
    seq.addBlock(mr.makeLabel('INC','REP', 1)); 
end
%% Timing check and report
outfile = [outfile 'linescan_tse_' dwType];
disp(outfile)

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

report=seq.testReport();

fprintf('Test Report listing follows:\n');
fprintf([report{:}]);
fprintf('\n');

[total_energy, peak_pwr, rf_rms] = seq.calcRfPower('windowDuration',min_seq_duration);
disp(['Sequence duration: ' sprintf('%0.3f',min_seq_duration)])
disp(['Total energy: ' sprintf('%0.3e',total_energy)])
disp(['Peak power: ' sprintf('%0.3e',peak_pwr/gamma^2) ' [mT^2]'])
disp(['RF rms: ' sprintf('%0.3e',rf_rms/gamma) ' [T]'])
%%
seq.setDefinition('FOV',[fov fov sliceThickness]);
seq.setDefinition('Name','ls_tse');
seq.setDefinition('SlicePositions',0);
seq.setDefinition('SliceThickness',sliceThickness);
seq.setDefinition('SliceGap',0);
seq.setDefinition('TR',TR);
% seq.setDefinition('kSpaceCenterLine', Npe/2+1);
%% Plot time diagram
seq.plot('timeRange', [0 2*TR],'stacked',1,'timeDisp','ms');
%% Evaluate label settings more specifically (is this necessary?)
lbls=seq.evalLabels('evolution','adc');
lbl_names=fieldnames(lbls);
figure; nexttile; hold on;
for n=1:length(lbl_names)
    plot(lbls.(lbl_names{n}),'linewidth',2);
end
grid on
legend(lbl_names(:));
title(['evolution of labels/counters/flags']);
xlabel('adc number');
%% k-space trajectory calculation and plotting
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
nexttile; hold on; plot(t_ktraj, ktraj'); xlim([0,min_seq_duration]);xlabel("t (s)") % plot the entire k-space trajectory
nexttile;plot(ktraj(1,:),ktraj(2,:),'b',...
    ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot

%% Check frequencies
FB(1).freq=590;
FB(1).bw=100;
FB(2).freq=1140;
FB(2).bw=220;
ax = nexttile;
[R, Rax, f] = seq.gradSpectrum(FB,[3000],false);
plot(f,Rax);
hold on; plot(f,R); % sos
xlabel('frequency (Hz)');
for i=1:length(FB)
    xline(FB(i).freq,'-');
    xline(FB(i).freq-FB(i).bw/2,'--');
    xline(FB(i).freq+FB(i).bw/2,'--');
end
disp(['Forbidden frequency power: ' sprintf('%0.3e', (sum(R(f<(FB(1).freq+FB(1).bw/2) & f>(FB(1).freq-FB(1).bw/2))) + sum(R(f<(FB(2).freq+FB(2).bw/2) & f>(FB(2).freq-FB(2).bw/2))))/(sum(FB(1).freq+FB(1).bw/2 & f>(FB(1).freq-FB(1).bw/2)) + sum(FB(2).freq+FB(2).bw/2 & f>(FB(2).freq-FB(2).bw/2))) )])
legend({'Gx','Gy','Gz','Gtot'});


%% Generate .seq file for MR scanner

% seq.write([outfile '.seq'])