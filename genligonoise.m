%% Generate simulated GW detector noise
% File containing target sensitivity curve (first column is frequency and
% second column is square root of PSD).
targetSens = load('iLIGOSensitivity.txt');
%%
% Plot the target sensitivity.
loglog(targetSens(:,1),targetSens(:,2));
xlabel('Frequency (Hz)');
ylabel('Strain Sensitivity (1/\sqrt{Hz})');

%%
% Select pass band.
fLow = 40;%Hz
fHigh = 1024;%Hz
indxfCut = targetSens(:,1)<= fLow | targetSens(:,1)>= fHigh;
targetSens(indxfCut,2)=0;
hold on;
loglog(targetSens(:,1),targetSens(:,2));
snapnow;
%%
% Sampling frequency of the data to be generated (should be less than half
% of the maximum frequency in target PSD.
fs = 4096;%Hz

%% Design filter
% B = fir2(N,F,A) designs an Nth order linear phase FIR digital filter with
% the frequency response specified by vectors F and A and returns the
% filter coefficients in length N+1 vector B.  The frequencies in F must be
% given in increasing order with 0.0 < F < 1.0 and 1.0 corresponding to
% half the sample rate.
%%
% FIR filter order
filtOrdr = 100;
%%
% We only include frequencies up to the Nyquist frequency (half of sampling
% frequency) when designing the filter.
indxfCut = targetSens(:,1)<=fs/2;
% Truncate the target PSD to Nyquist frequency.
targetSens = targetSens(indxfCut,:);

%%
% Add 0 frequency and corresponding PSD value
% as per the requirement of FIR2. Similarly add Nyquist frequency.
if targetSens(1,1) > 0 
    addZero = 1;
else
    addZero = 0;
end
if targetSens(end,1) < fs/2
    addNyq = 1;
else
    addNyq = 0;
end
%%
if addZero
    targetSens = [[0,0];targetSens];
end
if addNyq
    targetSens = [targetSens;[fs/2,0]];
end
%%
% Obtain filter coefficients. 
b = fir2(filtOrdr,targetSens(:,1)/(fs/2),targetSens(:,2));

%%
% Compare target and designed quantities
%Get the impulse response
impDataNSamples = 2048;
impSample = floor(impDataNSamples/2);
impVec = zeros(1,impDataNSamples);
impVec(impSample)=1;
impResp = fftfilt(b,impVec);
%%
%Get the transfer function
designTf = fft(impResp);
%%
%Plot the magnitude of the filter transfer function.
figure;
kNyq = floor(impDataNSamples/2)+1;
posFreq = (0:(kNyq-1))*(1/(impDataNSamples/fs));
plot(posFreq,abs(designTf(1:kNyq)));
hold on;
plot(targetSens(:,1),targetSens(:,2));
ylabel('TF magnitude');
legend('Designed','Target');
xlabel('Frequency (Hz)');
snapnow;

%% Generate noise
% Pass a white noise sequence through the designed filter.
nDataSamples = 16384;
inputNoise = randn(1,nDataSamples);
outputNoise = fftfilt(b,inputNoise);
figure;
plot((0:(nDataSamples-1))/fs,outputNoise);
snapnow;
%%
% Estimate PSD of simulated noise. *Note*: Scaling may be off because of
% (a) factors involved in discrete version of Wiener-Khinchin theorem, and
% (b) factors involved in how pwelch defines PSD. This can be easily
% corrected by multiplying with an overall factor (exercise: work out the
% factor!).
[pxx,f] = pwelch(outputNoise,2048,[],2048,fs);
figure;
plot(f,sqrt(pxx));
xlabel('Frequency (Hz)');
ylabel('[PSD]^{1/2}');
snapnow;