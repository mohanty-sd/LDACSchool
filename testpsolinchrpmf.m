%% Test script for PSO on the linear chirp matched filtering problem
% Generate data realization for signal and data parameters specified in the
% order below.
% Number of samples in data;
% Sampling frequency;
% SNR (signal amplitude); 
% f0;
% f1;
% initial phase;
% time of arrival;
% signal length.
nSamples = 2048;
fs = 2048;
sigLen = 0.25;
f0 = 384;%Hz
f1 = 39;%Hz^2
snr = 20;
[dataVec,timeVec]=genlinchrpwgndata(nSamples,fs,snr,f0,f1,pi/3.3,0.3,sigLen);
%% Spectrogram of data realization
[S,F,T]=spectrogram(dataVec,64,63,[],2048);
imagesc(T,F,abs(S)); axis xy;
snapnow;
%% Run PSO
% The fitness function called is LINCHIRPMFFITFUNC. The FFT of the data
% realization is passed through the fitness function input parameter
% structure (2nd input argument). Make sure that the search range includes
% the signal parameters as specified in the generation of the data
% realization. The first search parameter is f0 and the second is f1.
fftData = fft(dataVec);
ffparams = struct('rmin',[50,5],...
                     'rmax',[500,50],...
                     'fftData',fftData,...
                     'timeVec',timeVec,...
                     'fs',fs,...
                     'sigLen',sigLen);
%%
% Fitness function handle.
fitFuncHandle = @(x) linchirpmffitfunc(x,ffparams);
%%
% Call PSO.
psoOut = ldacpso(fitFuncHandle,2);

%% Estimated parameters
% Best standardized and real coordinates found.
stdCoord = psoOut.bestLocation;
[~,realCoord] = fitFuncHandle(stdCoord);
fprintf('Estimated f0=%f; Real f0=%f\n',realCoord(1),f0);
fprintf('Estimated f1=%f; Real f1=%f\n',realCoord(2),f1);