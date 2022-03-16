%% 

clearvars; 
% close all;


%% Definitions

% % Old (Proton):
% SBW = 120000;
% Omega0 = 297.223*10^6;
% T2 = 5E-3;


% Deuterium: Need min 6 ppm ~ 6*50 Hz = 300 Hz. 
% Peaks: Water @ 4.8 ppm, Glucose @ 3.8 ppm, Glu + Gln @ 2.4 ppm, Lac @ 1.4 ppm. So would need only 3.4 ppm, but peaks are broad --> Use 6 ppm.
% Max SBW that can be achieved for 1 TI and a 32x32 matrix: ~ 542 Hz (however 550 +- 50 Hz not allowed at scanner. So use 500 Hz)
% I calculated this like that: The gyro-ratio (1H/2H) = 26.752/4.107 ~ 6.51. So if we measure a 32x32 matrix, this corresponds to a 6.51*(32x32) ~ 210x210 matrix
% for protons. Now I chose that in poet with 1 TI, and got max 542 Hz SBW.
% SBW = 500;
% SBW = 10000; % Just for now to test

AcqTime = 1e-2;
% AcqTime = 1e1; % Just for now to test
AcqTimeVec = 4e-3:2e-3:5e-2;
% AcqTimeVec = 1e0;


Omega0 = 45.755 * 10^6;     % gamma_2H / 2 / pi * 7 T = 41.07 *1e6 / 2 / pi * 7T
T2 = 30E-3;
Freq_ppm = [1.5 2.5];
FreqsForAcqTimeSim_ppm = {[1.5 1.75],[1.5 2.0],[1.5 2.5]};
AmpsGroundTruth = [1 1];
SNR_Time_At500Hz = 5;
MonteCarloRealizations = 500;
SBWVec = 200:100:10000;
% SBWVec = 500:100:500;
SBW = 500;


%% Relative Error vs SBW


Index = 0;
RelativeMeanError = zeros([2 numel(SBWVec)]);
SNROfAmpEstimates = RelativeMeanError;
PlotDataSBWDep(numel(SBWVec))=struct('ppm',[],'Time',[],'FID',[],'Spec',[],'FitTot',[],'FitComps',[]);
for CurSBW = SBWVec
    Index = Index + 1;
    vecSize = round(AcqTime/(1/CurSBW));     % 10 ms / dt
    SNR_Time_Abs = SNR_Time_At500Hz / sqrt(CurSBW/500);

    [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(Freq_ppm,4.65,0,0,T2,AmpsGroundTruth,0,1/CurSBW,vecSize,Omega0);
 
    [MaxAmp,NoiseStd] = CalcNoise(FID_GroundTruth,SNR_Time_Abs);

    [RelativeMeanError(:,Index),SNROfAmpEstimates(:,Index),PlotDataSBWDep(Index)] = MonteCarloSimFitDeuteriumFID(MonteCarloRealizations,FID_GroundTruth,NoiseStd,T2,Omega0,Freq_ppm,Time,AmpsGroundTruth,PlotDataSBWDep(Index));
    PlotDataSBWDep(Index).ppm = ppm;
    
end
figure; plot(SBWVec,100*RelativeMeanError(1,:))
xlabel('SBW [Hz]'), ylabel('Relative Error [%]')

% Plot Data & Fit of one realization
PlotData(PlotDataSBWDep(4))

%% Relative Error vs AcqTime

RelativeMeanError = zeros([2 numel(AcqTimeVec) numel(FreqsForAcqTimeSim_ppm)]);
SNROfAmpEstimates = RelativeMeanError;
SNR_Time_Abs = SNR_Time_At500Hz / sqrt(SBW/500);
PlotDataAcqTimeDep(numel(AcqTimeVec),numel(FreqsForAcqTimeSim_ppm))=struct('ppm',[],'Time',[],'FID',[],'Spec',[],'FitTot',[],'FitComps',[]);


for FreqInd = 1:numel(FreqsForAcqTimeSim_ppm)
    CurFreq_ppm = FreqsForAcqTimeSim_ppm{FreqInd};
    Index = 0;
    for CurAcqTime = AcqTimeVec
        Index = Index + 1;
        vecSize = round(CurAcqTime/(1/SBW));     % 10 ms / dt

        [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(CurFreq_ppm,4.65,0,0,T2,AmpsGroundTruth,0,1/SBW,vecSize,Omega0);

        [MaxAmp,NoiseStd] = CalcNoise(FID_GroundTruth,SNR_Time_Abs);

        [RelativeMeanError(:,Index,FreqInd),SNROfAmpEstimates(:,Index,FreqInd),PlotDataAcqTimeDep(Index,FreqInd)] = MonteCarloSimFitDeuteriumFID(MonteCarloRealizations,FID_GroundTruth,NoiseStd,T2,Omega0,CurFreq_ppm,Time,AmpsGroundTruth,PlotDataAcqTimeDep(Index,FreqInd));
        PlotDataAcqTimeDep(Index,FreqInd).ppm = ppm;
    end
end
figure; plot(AcqTimeVec*1000,squeeze(100*RelativeMeanError(1,:,:)))
xlabel('Acquisition Duration [ms]'), ylabel('Relative Error [%]')
legend('0.25 ppm Peak Distance','0.5 ppm Peak Distance','1.0 ppm Peak Distance')

% Plot Data & Fit of one realization
PlotData(PlotDataAcqTimeDep(end))



%% Functions

function [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(Chemshift,DeltaFrequency,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq)

    Tmp = cell([1 numel(Chemshift)]);
    for ii = 1:numel(Chemshift)
        [Tmp{ii},ppm] = Simulate_FID_Spectra(Chemshift(ii),DeltaFrequency,phase0,AcqDelay,T2,S_0(ii),SNR,dwelltime,vecSize,LarmorFreq);   % Lac  % sqrt(2): We add both spectra, so the noise adds up
    end

    Time = Tmp{1}(1,:);
    ppm = ppm(1,:);

    FID_GroundTruth = [];
    for ii = 1:numel(Tmp)
       FID_GroundTruth = cat(3,FID_GroundTruth,Tmp{ii}(2,:)); 
    end
    FID_GroundTruth = squeeze(FID_GroundTruth);
    clear Tmp


end


function [MaxAmp,NoiseStd] = CalcNoise(FID_GroundTruth,SNR_Time_Abs)
    % Calc Noise
    MaxAmp = max(real(FID_GroundTruth(:)));
    NoiseStd = MaxAmp / (sqrt(2)*SNR_Time_Abs);      % SNR = MaxAmp / (sqrt(2)*std); --> std = MaxAmp / (sqrt(2)*SNR)
    NoiseStd = NoiseStd/sqrt(size(FID_GroundTruth,2)); % Because we add several spectral components --> Noise will add up!
end

function [RelativeMeanError,SNROfAmpEstimates,DataToPlot] = MonteCarloSimFitDeuteriumFID(MonteCarloRealizations,FID_GroundTruth,NoiseStd,T2,Omega0,Freq_ppm,Time,AmpsGroundTruth,DataToPlot)

    AmpError = zeros([MonteCarloRealizations size(FID_GroundTruth,2)]);
    Freq_ppm_tmp = transpose(Freq_ppm);
    T2_tmp = repmat(T2,[size(FID_GroundTruth,2) 1]);
    for ii = 1:MonteCarloRealizations
        Noise = NoiseStd * (randn(size(FID_GroundTruth)) +  1i*randn(size(FID_GroundTruth)));
        FID = FID_GroundTruth + Noise;

        FID_sum = transpose(squeeze(sum(FID,2)));


        % Fit FID with known frequencies and T2s
        y = @(T2Value,t,Omegas,x) sum(x.*exp(-t./T2Value).*exp(-1i*Omegas.*t),1);


        RealOmega = Omega0 * (1 + (Freq_ppm_tmp - 4.65)/10^6) * 2*pi;
        RealOmega = RealOmega - Omega0 *2*pi;

        OLS = @(x) sum(abs((y(T2_tmp,Time,RealOmega,x) - FID_sum)).^2);           % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',500000, 'MaxIter',9000000);
        A = fminsearch(OLS, zeros([size(FID_GroundTruth,2) 1]), opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function
        A = transpose(A);

        AmpError(ii,:) = A - AmpsGroundTruth;

    end

    DataToPlot.Time = Time;
    DataToPlot.FID = FID;
    DataToPlot.Spec = fftshift(fft(FID,[],1),1);
    DataToPlot.FitTot = y(repmat(T2,[size(FID_GroundTruth,2) 1]),Time,RealOmega,transpose(A));                            % Calculate function with estimated parameters

    for ii = 1:size(FID_GroundTruth,2)
        Zerr = zeros([size(FID_GroundTruth,2) 1]);
        Zerr(ii) = A(ii);
        DataToPlot.FitComps{ii} = y(repmat(T2,[size(FID_GroundTruth,2) 1]),Time,RealOmega,Zerr);
    end


    % Relative Error & SNR
%     RelativeError = AmpError ./ AmpsGroundTruth;
    RelativeMeanError = mean(abs(AmpError ./ AmpsGroundTruth));

    SNROfAmpEstimates = AmpsGroundTruth ./ std(AmpError);


end


function PlotData(DataToPlot)
    LegendTitle = {'Data','Fit'};
    

    figure; 

    
    % Plot Individual Spectra
    subplot(2,2,1)
    hold on;
    for ii = 1:numel(DataToPlot.FitComps)
        Linee = plot(DataToPlot.ppm,real(fftshift(fft(DataToPlot.FitComps{ii}))));
        colour = get(Linee,'Color');
        scatter(transpose(DataToPlot.ppm),real(DataToPlot.Spec(:,ii)),40,colour,'filled')
    end
    hold off;
    xlabel('Chemical Shift [ppm]'), ylabel('Signal [a.u.]')
    legend(LegendTitle)
    title('Separate Specs + Fits')
    
    % Plot Individual FIDs
    % Scatter Points
    subplot(2,2,2) 
    hold on
    
    
    for ii = 1:numel(DataToPlot.FitComps)
        Linee = plot(1000*DataToPlot.Time,real(DataToPlot.FitComps{ii}));
        colour = get(Linee,'Color');        
        scatter(1000*DataToPlot.Time,real(DataToPlot.FID(:,ii)),60,colour,'filled')
    end
    hold off
    % Label
    xlabel('Time [ms]'), ylabel('Signal [a.u.]')
    legend(LegendTitle)
    title('Separate FIDs + Fits')
    
    
    % Plot Sum Spectrum + Fit
    subplot(2,2,3); 
    scatter(DataToPlot.ppm,real(sum(DataToPlot.Spec,2)),60,'b','filled')
    hold on
    % Line
    plot(DataToPlot.ppm,real(fftshift(fft(DataToPlot.FitTot))),'b')
    hold off
    xlabel('Chemical Shift [ppm]'), ylabel('Signal [a.u.]')
    legend('Lac+Glx','Fit')
    title('Sum Spec + Fit')
    
    
    % Plot Sum FID + Fit
    subplot(2,2,4); 
    scatter(1000*DataToPlot.Time(1,:),real(sum(DataToPlot.FID,2)),60,'b','filled')
    hold on
    % Line
    plot(1000*DataToPlot.Time(1,:),real(DataToPlot.FitTot),'b')
    xlabel('Time [ms]'), ylabel('Signal [a.u.]')
    hold off
    legend('Lac+Glx','Fit')
    title('Sum FID + Fit')

end
