%% Calculate the TDOA measurment from the two spectral signals
function [TDOA,GPHATi,E] = TDOAmeasurment(b1, b2, Para)

b1f = fft(b1,Para.NFFT); b1f = b1f(1:Para.NFFT/2+1);
b2f = fft(b2,Para.NFFT); b2f = b2f(1:Para.NFFT/2+1);

% voice band energy calculation
temp = b1f(Para.voiceFre(1)+1:Para.voiceFre(2)+1);
E = temp'*temp;
temp = b2f(Para.voiceFre(1)+1:Para.voiceFre(2)+1);
E = E+temp'*temp;

% apply PHAT-GCC
b2fc=conj(b2f);
neuma=(b1f).*(b2fc);
deno=abs((b1f).*(b2fc))+eps;
GPHAT=neuma./deno;
% GPHAT2 = [zeros(1);GPHAT;zeros(1)];

tau_array = Para.tau_array;
expwtau = Para.expwtau;
Ntau = Para.Ntau;

temp = repmat(GPHAT,[1 Ntau]).*expwtau;
GPHATi = real(sum(temp)); GPHATi = GPHATi(:);

% % smooth the spectrum
% GPHATi = 0.8*GPHATi+0.2*GPHATi_previous;

% Now find the peaks from the GCC functions over candidates taus, and save
% them as TDOA measurments.
NdoaMax = Para.NdoaMax; % find at most NdoaMax measurements at each frame t

radiusN = Para.radiusN;
elem_x = Para.elem_x;

Threshold = Para.Threshold;

i = 1;
ind_array = [];
maxV_array = [];

% colorarray = 'rbgcmk';
% figure(100);
% clf
% hold on;
while (i<=NdoaMax)
%     plot(GPHATi,colorarray(i*2-1));
    % find the maximum value
    [maxV,ind] = max(GPHATi);
    affectIndex = ind+(-radiusN:radiusN);
    useIndex = affectIndex>1 & affectIndex<Ntau;
    realIndex = affectIndex(useIndex);
    GPHATi(realIndex) = GPHATi(realIndex)-maxV*elem_x(useIndex);
%     plot(GPHATi,colorarray(i*2));
    if i==1, % find at least one TDOA measurement
        ind_array = [ind_array;ind];
        maxV_array = [maxV_array; maxV];
    else
        if maxV>Threshold && maxV>0.5*maxV_array(1),
            ind_array = [ind_array;ind];
            maxV_array = [maxV_array; maxV];
        else
            break
        end
    end
    i=i+1;
end
hold off

TDOA = [-tau_array(ind_array)' maxV_array];



