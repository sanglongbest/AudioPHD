function scores = gccscores(a,settings)
%GCCSCORES -
%
%    This function...
%
%    scores = GCCSCORES(a,settings)
%
%    Input:
%    a - matrix
%    settings - struct that must have...
%
%    Output:
%    scores - two-dimensional cell-array containg...

%Calculating matrices containing the STFT of all frames:
fftframes = cell(1,settings.mm);
for k = settings.channels
    [fftframes{k}, f, t] = stft(a(k,:), settings.frameSize, settings.dx, settings.frameSize*2, settings.sr);
    %fftframes{k} = fftframes{k}(1:size(fftframes{k},1)/2+1,:);
end
clear frames;

% apply PHAT-GCC
fftframes{2}=conj(fftframes{2});
neuma=(fftframes{1}).*(fftframes{2});
deno=abs((fftframes{1}).*(fftframes{2}))+eps;
GPHAT=neuma./deno;
% GPHAT2 = [zeros(1);GPHAT;zeros(1)];

tau_array = [-1:0.1:1];
expwtau = exp(sqrt(-1)*(1./f)'*tau_array);
Ntau = size(tau_array,2)

temp = repmat(GPHAT,[1 Ntau]).*expwtau;
GPHATi = real(sum(temp)); GPHATi = GPHATi(:);
scores = GPHATi








% %Getting indeces indu only for the pairs that we have to calculate:
% ii = zeros(1,settings.mm);
% ii(settings.channels) = 1;
% ind = logical(ii'*ii);
% indu = triu(ind);
% 
% %Cross-correlation for all pairs with indeces indu:
% scores = cell(settings.mm);
% for k = find(indu)'
%     [ii,jj] = ind2sub(size(ind),k);
%     tmp = fftframes{jj}.*conj(fftframes{ii});
%     tmp = fftshift(ifft(tmp.*settings.wf(tmp)),1);
%     scores{ii,jj} = tmp(round(end/2)-settings.sw:round(end/2)+settings.sw,:);
% end
% clear fftframes;
% 
% %Filling in the pairs not in indu by mirroring the data:
% tmp = cellfun(@flipud,scores','UniformOutput',false);
% scores(tril(ind,-1)) = tmp(tril(ind,-1));
end