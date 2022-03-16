%% 

clearvars; 
% close all;


%% Definitions

% % Old (Proton):
% SBW = 120000;
% Omega0 = 297.223*10^6;
% T2 = 5E-3;


% Deuterium: Need min 6 ppm ~ 6*50 Hz = 300 Hz. 
% Max SBW that can be achieved for 1 TI and a 32x32 matrix: ~ 542 Hz (however 550 +- 50 Hz not allowed at scanner. So use 500 Hz)
% I calculated this like that: The gyro-ratio (1H/2H) = 26.752/4.107 ~ 6.51. So if we measure a 32x32 matrix, this corresponds to a 6.51*(32x32) ~ 210x210 matrix
% for protons. Now I chose that in poet with 1 TI, and got max 542 Hz SBW.
% SBW = 500;
% SBW = 10000; % Just for now to test

AcqTime = 1e-2;
% AcqTime = 1e1; % Just for now to test
AcqTimeVec = 4e-3:2e-3:2e-1;
% AcqTimeVec = 1e0;


Omega0 = 45.755 * 10^6;     % gamma_2H / 2 / pi * 7 T = 41.07 *1e6 / 2 / pi * 7T
T2 = 30E-3;
Freq_ppm = [1.5 2.5];
AmpsGroundTruth = [1 1];
SNR_Time_At500Hz = 5;
MonteCarloRealizations = 1000;
SBWVec = 200:100:10000;
% SBWVec = 500:100:500;
SBW = 500;


%% Relative Error vs SBW


% Index = 0;
% RelativeMeanError = zeros([2 numel(SBWVec)]);
% SNROfAmpEstimates = RelativeMeanError;
% PlotDataSBWDep(numel(SBWVec))=struct('ppm',[],'Time',[],'FID',[],'Spec',[],'FitTot',[],'FitComp1',[],'FitComp2',[]);
% for CurSBW = SBWVec
%     Index = Index + 1;
%     vecSize = round(AcqTime/(1/CurSBW));     % 10 ms / dt
%     SNR_Time_Abs = SNR_Time_At500Hz / sqrt(CurSBW/500);
% 
%     [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(Freq_ppm,4.65,0,0,T2,AmpsGroundTruth,0,1/CurSBW,vecSize,Omega0);
%  
%     [MaxAmp,NoiseStd] = CalcNoise(FID_GroundTruth,SNR_Time_Abs);
% 
%     [RelativeMeanError(:,Index),SNROfAmpEstimates(:,Index),PlotDataSBWDep(Index)] = MonteCarloSimFitDeuteriumFID(MonteCarloRealizations,FID_GroundTruth,NoiseStd,T2,Omega0,Freq_ppm,Time,AmpsGroundTruth,PlotDataSBWDep(Index));
%     PlotDataSBWDep(Index).ppm = ppm;
%     
% end
% figure; plot(SBWVec,RelativeMeanError(1,:))


%% Relative Error vs AcqTime

Index = 0;
RelativeMeanError = zeros([2 numel(AcqTimeVec)]);
SNROfAmpEstimates = RelativeMeanError;
SNR_Time_Abs = SNR_Time_At500Hz / sqrt(SBW/500);
PlotDataAcqTimeDep(numel(SBWVec))=struct('ppm',[],'Time',[],'FID',[],'Spec',[],'FitTot',[],'FitComp1',[],'FitComp2',[]);
for CurAcqTime = AcqTimeVec
    Index = Index + 1;
    vecSize = round(CurAcqTime/(1/SBW));     % 10 ms / dt

    [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(Freq_ppm,4.65,0,0,T2,AmpsGroundTruth,0,1/SBW,vecSize,Omega0);
 
    [MaxAmp,NoiseStd] = CalcNoise(FID_GroundTruth,SNR_Time_Abs);

    [RelativeMeanError(:,Index),SNROfAmpEstimates(:,Index),PlotDataAcqTimeDep(Index)] = MonteCarloSimFitDeuteriumFID(MonteCarloRealizations,FID_GroundTruth,NoiseStd,T2,Omega0,Freq_ppm,Time,AmpsGroundTruth,PlotDataAcqTimeDep(Index));
    
end
figure; plot(AcqTimeVec*1000,RelativeMeanError(1,:))






%% Functions

function [FID_GroundTruth,Time,ppm] = SimulateDeuteriumFIDs(Chemshift,DeltaFrequency,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq)
    [Tmp{1},ppm] = Simulate_FID_Spectra(Chemshift(1),DeltaFrequency,phase0,AcqDelay,T2,S_0(1),SNR,dwelltime,vecSize,LarmorFreq);   % Lac  % sqrt(2): We add both spectra, so the noise adds up
    [Tmp{2},ppm] = Simulate_FID_Spectra(Chemshift(2),DeltaFrequency,phase0,AcqDelay,T2,S_0(2),SNR,dwelltime,vecSize,LarmorFreq);   % Glx 

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

    AmpError = zeros([MonteCarloRealizations 2]);
    for ii = 1:MonteCarloRealizations
        Noise = NoiseStd * (randn(size(FID_GroundTruth)) +  1i*randn(size(FID_GroundTruth)));
        FID = FID_GroundTruth + Noise;

        FID_sum = transpose(squeeze(sum(FID,2)));


        % Fit FID with known frequencies and T2s
        y = @(T2Value,t,Omegas,x) x(1)*exp(-t/T2Value(1)).*exp(-1i*Omegas(1)*t) + x(2)*exp(-t/T2Value(2)).*exp(-1i*Omegas(2)*t);

        RealOmega = Omega0 * (1 + (Freq_ppm - 4.65)/10^6) * 2*pi;
        RealOmega = RealOmega - Omega0 *2*pi;

        OLS = @(x) sum(abs((y([T2,T2],Time,RealOmega,x) - FID_sum)).^2);           % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',500000, 'MaxIter',9000000);
        A = fminsearch(OLS, [0 0], opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function

        AmpError(ii,:) = A - AmpsGroundTruth;

    end

    DataToPlot.Time = Time;
    DataToPlot.FID = FID;
    DataToPlot.Spec = fftshift(fft(FID,[],1),1);
    DataToPlot.FitTot = y([T2,T2],Time,RealOmega,A);                            % Calculate function with estimated parameters

    DataToPlot.FitComp1 = y([T2,T2],Time,RealOmega,[A(1) 0]);
    DataToPlot.FitComp2 = y([T2,T2],Time,RealOmega,[0 A(2)]);



    % Relative Error & SNR
%     RelativeError = AmpError ./ AmpsGroundTruth;
    RelativeMeanError = mean(abs(AmpError ./ AmpsGroundTruth));

    SNROfAmpEstimates = AmpsGroundTruth ./ std(AmpError);


end


function PlotData(DataToPlot)
    LegendTitle = {'Lac','LacFit','Glx','GlxFit'};
    

    figure; 

    
    % Plot Individual Spectra
    subplot(2,2,1)
    scatter(transpose(DataToPlot.ppm),real(DataToPlot.Spec(:,1)),'b','filled'), 
    hold on;
    plot(DataToPlot.ppm,real(fftshift(fft(DataToPlot.FitComp1))),'b')
    scatter(transpose(DataToPlot.ppm),real(DataToPlot.Spec(:,2)),'r','filled'), 
    plot(DataToPlot.ppm,real(fftshift(fft(DataToPlot.FitComp2))),'r')
    hold off;
    xlabel('Chemical Shift [ppm]'), ylabel('Signal [a.u.]')
    legend(LegendTitle)
    title('Separate Specs + Fits')
    
    % Plot Individual FIDs
    % Scatter Points
    subplot(2,2,2) 
    scatter(1000*DataToPlot.Time(1,:),real(DataToPlot.FID(:,1)),60,'b','filled')
    hold on
    % Line
    plot(1000*DataToPlot.Time(1,:),real(DataToPlot.FitComp1),'b')
    scatter(1000*DataToPlot.Time(1,:),real(DataToPlot.FID(:,2)),60,'r','filled')
    plot(1000*DataToPlot.Time(1,:),real(DataToPlot.FitComp2),'r')
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
