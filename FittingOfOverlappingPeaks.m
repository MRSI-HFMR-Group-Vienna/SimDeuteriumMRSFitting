%% 

clearvars; 
% close all;


%%

% % Old (Proton):
% SBW = 120000;
% Omega0 = 297.223*10^6;
% T2 = 5E-3;


% Deuterium: Need min 6 ppm ~ 6*50 Hz = 300 Hz. 
% Max SBW that can be achieved for 1 TI and a 32x32 matrix: ~ 542 Hz (however 550 +- 50 Hz not allowed at scanner. So use 500 Hz)
% I calculated this like that: The gyro-ratio (1H/2H) = 26.752/4.107 ~ 6.51. So if we measure a 32x32 matrix, this corresponds to a 6.51*(32x32) ~ 210x210 matrix
% for protons. Now I chose that in poet with 1 TI, and got max 542 Hz SBW.
SBW = 500;
SBW = 2000; % Just for now to test

Omega0 = 45.755 * 10^6;     % gamma_2H / 2 / pi * 7 T = 41.07 *1e6 / 2 / pi * 7T
T2 = 30E-3;
vecSize = round(1e-2/(1/SBW));     % 10 ms / dt
Freq_ppm = [1.5 2.5];
Amps = [1 1];
SNR = 10;

[Tmp{1},ppm] = Simulate_FID_Spectra(Freq_ppm(1),4.65,0,0,T2,Amps(1),SNR*sqrt(2),1/SBW,vecSize,Omega0);   % Lac  % sqrt(2): We add both spectra, so the noise adds up
[Tmp{2},ppm] = Simulate_FID_Spectra(Freq_ppm(2),4.65,0,0,T2,Amps(2),SNR*sqrt(2),1/SBW,vecSize,Omega0);   % Glx 

Time = Tmp{1}(1,:);
ppm = ppm(1,:);

FID = [];
for ii = 1:numel(Tmp)
   FID = cat(3,FID,Tmp{ii}(2,:)); 
end
FID_sum = squeeze(sum(FID,3));
FID = squeeze(FID);
clear Tmp

Spec = fftshift(fft(FID,[],1),1); 




%% Fit FID with known frequencies and T2s


y = @(T2Value,t,Omegas,x) x(1)*exp(-t/T2Value(1)).*exp(-1i*Omegas(1)*t) + x(2)*exp(-t/T2Value(2)).*exp(-1i*Omegas(2)*t);

RealOmega = Omega0 * (1 + (Freq_ppm - 4.65)/10^6) * 2*pi;
RealOmega = RealOmega - Omega0 *2*pi;

OLS = @(x) sum(abs((y([T2,T2],Time,RealOmega,x) - FID_sum)).^2);           % Ordinary Least Squares cost function
opts = optimset('MaxFunEvals',500000, 'MaxIter',9000000);
A = fminsearch(OLS, [1 1], opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function
fcnfit = y([T2,T2],Time,RealOmega,A);                            % Calculate function with estimated parameters

Error = A - Amps;
% figure; scatter(Time,real(FID_sum)); hold on; plot(Time,real(fcnfit)); hold off;




%% Plot

LegendTitle = {'Lac','Glx'};


% Plot Individual Spectra
figure; 
plot(ppm,real(Spec)), xlabel('Chemical Shift [ppm]'), ylabel('Signal [a.u.]'), title('Spec'), legend(LegendTitle)


% Plot Individual FIDs
% Scatter Points
figure; scatter(1000*Time(1,:),real(FID(:,1)),60,'b','filled')
hold on
% Line
plot(1000*Time(1,:),real(FID(:,1)),'b')
scatter(1000*Time(1,:),real(FID(:,2)),60,'r','filled')
plot(1000*Time(1,:),real(FID(:,2)),'r')
hold off
% Label
xlabel('Time [ms]'), ylabel('Signal [a.u.]')


% Plot Sum Spectrum + Fit
figure; scatter(ppm,real(sum(Spec,2)),60,'b','filled')
hold on
% Line
plot(ppm,real(fftshift(fft(fcnfit))),'r')
hold off


% Plot Sum FID + Fit
figure; scatter(1000*Time(1,:),real(FID_sum),60,'b','filled')
hold on
% Line
plot(1000*Time(1,:),real(fcnfit),'r')
xlabel('Time [ms]'), ylabel('Signal [a.u.]')
hold off
