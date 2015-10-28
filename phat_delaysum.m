function delays = phat_delaysum(waves, fs, refIdx, margin)
%y = phat_delaysum(x, fs, refIdx, vMargin) returns a waveform obtained by
%delaying and summing the supplied input waves. Cross-channel time-delay 
%estimation (TDE) is performed using PHAT-weighted generalized 
%cross-correlation.
%
%Arguments:
%   x: cell array containing the input waveforms.
%   fs: sampling rate (in sps).
%   refIdx: index of the reference channel.
%   margin: search range for the maximum correlation lag (in seconds).

%This script is based on ArrayProcessor.m written by Tuomo Pirinen in spring
% 2004. It was modified by Marc Ferras in June 2005 to implement PHAT-GCC
%time-delay estimation.

%Switch waves so that the reference is always at index 1
tmp=waves{1};
waves{1}=waves{refIdx};
waves{refIdx}=tmp;

nChannels = length(waves) ; %Number of input channels
nSamples=length(waves{1});

marginSamples = round(margin*fs);  %Margin in samples

output = zeros(nSamples,1) ; %Initialize output vector
delays=zeros(1,nChannels); %Initizalize TDE vector

%Find 2^x such that 2^x>2*nSamples
fftSize=2;
while fftSize<2*nSamples
        fftSize=fftSize*2;
end
fftSize2=floor(fftSize/2);

for(k=1:nChannels)
  %GCC-PHAT
  fft1=fft(waves{k},fftSize);
  fft2=fft(waves{1},fftSize);
  G12=fft1.*conj(fft2);
  denom=max(abs(G12),1e-6);
  G=G12./denom;
  f=real(ifft(G));
  %plot(f,'-b');
  g=fftshift(f);
  %hold on;
  %plot(g,'-r');
  
  [maxVal, maxIdx] = max(g(fftSize2+1-marginSamples:fftSize2+1+marginSamples));
  delays(k) = maxIdx - marginSamples - 1;  
  
%   mtxMatchSamples(k,1) = [ delays(k) + 1 ] ;
%   if(mtxMatchSamples(k,1)<1)
%     disp('Pre-padding')
%     padLen = 1 - delays(k) ;
%     waves{k} = [ zeros(padLen,1) ; waves{k} ] ;
%     mtxMatchSamples(k,1) = 1 ;
%     mtxMatchSamples(k,2) = mtxMatchSamples(k,2) + padLen ;
%   end
% 
%   mtxMatchSamples(k,2) =  mtxMatchSamples(k,1) + nSamples - 1 ;
%   if( mtxMatchSamples(k,2) > nSamples )
%     disp('Post-padding') ;
%     padLen = mtxMatchSamples(k,2) - nSamples + 1; 
%     waves{k} = [ waves{k} ; zeros(padLen,1)  ] ;    
%   end
%   
%   alignedWave = waves{k}(mtxMatchSamples(k,1):mtxMatchSamples(k,2));
%   output = output + alignedWave ;
end

delays;
%output = 0.8*output./max(abs(output)) ;

return