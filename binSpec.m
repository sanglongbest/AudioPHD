function [L,R, W, HopSize] = binSpec(lr, NFFT, Window_or_WindowSize, HopSize_or_HopPercent)

% [L,R, W, HopSize] = binSpec(lr, NFFT, Window_or_WindowSize,
% HopSize_or_HopPercent)
%
% Make a spectrogram of the left and right channels of the waveform
% LR.  LR(1,:) is the left channel, LR(2,:) is the right channel.
% NFFT is the size of the window to use and the number of samples
% in each window, HOP is the number of samples between window
% starts.  Note that this also chops off the DC and nyquist
% components.

if nargin<2, NFFT = 1024; end
if nargin<3, Window_or_WindowSize = NFFT;end
if nargin<4, HopSize_or_HopPercent=.25; end

if length(Window_or_WindowSize)==1,
 	Window=hamming(Window_or_WindowSize);
else
 	Window=Window_or_WindowSize(:);
end
W = length(Window);

if HopSize_or_HopPercent>1,
	HopSize = round(HopSize_or_HopPercent);
else
	HopSize = fix(W.*HopSize_or_HopPercent);
end


if(ndims(lr) == 3)
  % Already specgrammed with DC and Nyquist Rate deleted
  L = lr(:,:,1);
  R = lr(:,:,2);
else
  L = stft(lr(1,:), NFFT, Window_or_WindowSize, HopSize_or_HopPercent);
  R = stft(lr(2,:), NFFT, Window_or_WindowSize, HopSize_or_HopPercent);

  L = L(2:end-1, :);
  R = R(2:end-1, :);
end
