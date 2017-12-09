function [dataVec, timeVec] = genlinchrpwgndata(nSamples, fs, ...
                                                sigAmplitude, f0, f1,...
                                                initialPhase, timeOfArrival,...
                                                sigLen)
%Generate white Gaussian noise plus linear chirp data realization
%[Y,T] = GENLINCHIRPWGNDATA(N,Fs,A,F0,F1,P0,TOA,SL)
% Generate a data realization Y containing a linear chirp signal of length
% SL, frequency parameters F0 (Hz), F1 (Hz^2), initial phase P0, time of arrival (in
% seconds) TOA, and amplitude A added a realization of i.i.d. Normal random
% variable sequence (White Gaussian Noise) with zero mean and unit
% variance. The sampling times are returned in vector T.

%Soumya D. Mohanty, Dec 2017

%Sampling times
timeVec = (0:(nSamples-1))/fs;

%% 
%Length of data in seconds
dataLen = timeVec(end);

%%
% Find the samples where the signal should start and end. Store the
% corresponding time samples.
startSample = floor(timeOfArrival*fs);
endSample = floor((timeOfArrival+sigLen)*fs);
sigTimeVec = timeVec(startSample:endSample);
%% Generate the signal
%First allocate an empty array that is as long as the data. Then load the
%signal between start and end samples corresponding to the time of arrival
%and the end of the signal.
sigVec = zeros(1,nSamples);
sigVec(startSample:endSample) = sin(2*pi*(f0*(sigTimeVec-timeOfArrival)+...
                   f1*(sigTimeVec-timeOfArrival).^2)+...
                   initialPhase);
               
% Normalize the signal and multiply with the amplitude parameter.
sigVec = sigVec/norm(sigVec);
sigVec = sigAmplitude*sigVec;

% Noise realization
noiseVec = randn(1,nSamples);

% Add signal and noise
dataVec = sigVec+noiseVec;


