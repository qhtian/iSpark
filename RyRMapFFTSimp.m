function [TTPower,I_fft_amp_norm_smoothn]=RyRMapFFTSimp(I,xyt_dim,livH)
    s=2.0; % Original=2; %% Smooth s for fft amplitude.
    radius=3/xyt_dim(1);
    siz=size(I);
    
    %% 2D FFT amplitude.
    I=double(I);
    I_fft_amp=abs(fftshift(fft2(I)));
    
    %% Calculate its surrounding background and normalize amplitude.
    radius=floor(radius/2)*2+1;
    h=fspecial('disk',radius)>0;
    h=h-ordfilt2(h,1,[0,1,0;1,1,1;0,1,0],'zeros')==1;
    
    nthMedian=h>0;
    nthMedian=floor(sum(nthMedian(:))/2);
    
    Bgr=ordfilt2(I_fft_amp,nthMedian,h,'symmetric');
    
    I_fft_amp_norm=I_fft_amp./Bgr;
    I_fft_amp_norm_smoothn=smoothn(I_fft_amp_norm,s);
    clear('Bgr','nthMedian','h','radius','I_fft_amp')
    
    %% Get the center point
    I=I/max(I(:))*5;    
    % if isempty(x)
    x=findFFTPeak(I_fft_amp_norm_smoothn);
    if isempty(x)
        [y,x] = getpts(gcf);
    else
        y=x(1,2);x=x(1,1);
    end
    % else
    %     y=x(1,2);x=x(1,1);
    % end
    Centroid=cat(1,[x,y],[size(I,1)/2,size(I,2)/2]);
    dist2Center=sqrt((Centroid(2,2)-Centroid(1,2))^2+(Centroid(2,1)-Centroid(1,1))^2);
    
    %% Refine the centroid.
    cx1b=round(x(1)-dist2Center*0.2); if cx1b<1; cx1b=1; end
    cy1b=round(y(1)-dist2Center*0.2); if cy1b<1; cy1b=1; end
    cx1e=round(x(1)+dist2Center*0.2); if cx1e>siz(1); cx1e=siz(1); end
    cy1e=round(y(1)+dist2Center*0.2); if cy1e>siz(2); cy1e=siz(2); end
    
    centralArea=I_fft_amp_norm_smoothn(cx1b:cx1e,cy1b:cy1e);
    TTPower=max(centralArea(:));
    BW=centralArea==TTPower;
    x=find(sum(BW,2)==1);
    y=find(sum(BW,1)==1);
    x=x+cx1b-1; y=y+cy1b-1;
    
    %% Plotting
    if ishandle(livH) || livH
        if islogical(livH)
            livH=subplot(2,1,1);
            set(get(livH,'Parent'),'unit','normalized','Position',[0.15,0.1,0.75,0.8]);
        end
        
        image(I_fft_amp_norm_smoothn,'Parent',livH,'CDataMapping','scaled');
        axis(livH,'image','off');
        set(livH,'CLim',[0 5]);
        colormap(livH,CoolNature);
        hold(livH,'all');
        plot(livH,y,x,'*','Color',[0,0,0]);
        hold(livH,'off');
        title(livH,['FFT power: ' num2str(TTPower)])
        colorbar('peer',livH);

        livH=subplot(2,1,2);
        image(I,'Parent',livH,'CDataMapping','scaled');
        axis(livH,'image','off');
        colorbar('peer',livH);
        colormap(livH,CoolNature);
        try
           [~,hmax]=getRightImRange(I);
           set(livH,'CLim',[0,hmax]);
        catch
        end
        
        drawnow expose;
    end
    %% Final outputs.
    % x=[x,y];
end



function BlackToRed=CoolNature
    BlackToRed=[0 0 0;0 0 0.0666666701436043;0 0 0.133333340287209;0 0 0.200000002980232;...
        0 0 0.266666680574417;0 0 0.333333343267441;0 0 0.400000005960464;0 0 0.466666668653488;...
        0 0 0.533333361148834;0 0 0.600000023841858;0 0 0.666666686534882;0 0 0.733333349227905;...
        0 0 0.800000011920929;0 0 0.866666674613953;0 0 0.933333337306976;0 0 1;0 0.0625 1;...
        0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;...
        0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;...
        0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;...
        0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;...
        0.9375 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;...
        1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;...
        1 0.0625 0;1 0 0];
end

function [hmin,hmax]=getRightImRange(I)
    AUTO_THRESHOLD = 5000;
    pixcount=numel(I);
    limit = pixcount/10;
    threshold = pixcount/AUTO_THRESHOLD;
    nBins = 256;
    I=double(I);
    [histA,values]=hist(I(:),nBins);
    
    i = 0;
    found = false;
    while((~found)&&(i<numel(histA)))
        i=i+1;
        counts = histA(i);
        if (counts > limit)
            counts = 0;
        end
        found = counts > threshold;
    end
    hmin = values(i);
    
    i = numel(histA)+1;
    found = false;
    while((~found)&&(i>1))
        i=i-1;
        counts = histA(i);
        if (counts > limit)
            counts = 0;
        end
        found = counts > threshold;
    end
    hmax = values(i);
end

function [z,s,exitflag,Wtot] = smoothn(Raw,s,W)
    %% Test & prepare the variables
    Raw = single(Raw);
    sizeRaw = size(Raw);
    numOfElements = prod(sizeRaw); % number of elements
    if numOfElements<2, z = Raw; return, end
    
    % Smoothness parameter and weights
    if nargin==1; s=[]; W=ones(sizeRaw,'single'); elseif nargin==2; W=ones(sizeRaw,'single'); end
    if (~isequal(size(W),sizeRaw)) || (~isempty(s) && (~isscalar(s) || s<0)); z=Raw; return; end
    MaxIter = 100;      % "Maximal number of iterations" criterion
    TolZ = 1e-3;        % default value for TolZ
    
    % Weights. Zero weights are assigned to not finite values (Inf or NaN),
    IsFinite = isfinite(Raw);
    numOfFinites = nnz(IsFinite); % number of finite elements
    W = W.*IsFinite;        W = W/max(W(:));
    isweighted = any(W(:)<1);
    
    % Automatic smoothing?
    isauto = isempty(s);
    
    %% Creation of the Lambda tensor
    d = ndims(Raw);
    Lambda = zeros(sizeRaw,'single');
    for i = 1:d
        siz0 = ones(1,d);
        siz0(i) = sizeRaw(i);
        Lambda = bsxfun(@plus,Lambda,cos(pi*(reshape(1:sizeRaw(i),siz0)-1)/sizeRaw(i)));
    end
    Lambda = -2*(d-Lambda);
    if ~isauto, Gamma = 1./(1+s*Lambda.^2); end
    
    %% Upper and lower bound for the smoothness parameter
    N = sum(sizeRaw~=1); % tensor rank of the y-array
    hMin = 1e-6; hMax = 0.99;
    sMinBnd = (((1+sqrt(1+8*hMax.^(2/N)))/4./hMax.^(2/N)).^2-1)/16;
    sMaxBnd = (((1+sqrt(1+8*hMin.^(2/N)))/4./hMin.^(2/N)).^2-1)/16;
    
    %% Initialize before iterating
    Wtot = W;
    %--- Initial conditions for z
    if isweighted
        z = InitialGuess(Raw,IsFinite);
    else
        z = zeros(sizeRaw,'single');
    end
    
    z0 = z;
    Raw(~IsFinite) = 0; % arbitrary values for missing y-data
    
    tol = 1;
    numOfCurrIter = 0;
    %--- Error on p. Smoothness parameter s = 10^p
    errp = 0.1;    opt = optimset('TolX',errp);
    %--- Relaxation factor RF: to speedup convergence
    RF = 1 + 0.75*isweighted;
    
    %% Main iterative process
    
    %--- "amount" of weights (see the function GCVscore)
    aow = sum(Wtot(:))/numOfElements;           % 0 < aow <= 1
    while tol>TolZ && numOfCurrIter<MaxIter
        numOfCurrIter = numOfCurrIter+1;
        DCTy = dctn(Wtot.*(Raw-z)+z);
        if isauto && ~rem(log2(numOfCurrIter),1)
            fminbnd(@gcv,log10(sMinBnd),log10(sMaxBnd),opt);
        end
        z = RF*idctn(Gamma.*DCTy) + (1-RF)*z;
        
        % if no weighted/missing data => tol=0 (no iteration)
        tol = isweighted*norm(z0(:)-z(:))/norm(z(:));
        
        z0 = z; % re-initialization
    end
    exitflag = numOfCurrIter<MaxIter;
    
    
    %% Warning messages
    if nargout<3 && ~exitflag
        warning('MATLAB:smoothn:MaxIter',...
            ['Maximum number of iterations (' int2str(MaxIter) ') has ',...
            'been exceeded. Increase MaxIter option or decrease TolZ value.'])
    end
    
    %% GCV score
    %---
    function GCVscore = gcv(p)
        % Search the smoothing parameter s that minimizes the GCV score
        %---
        s = 10^p;
        Gamma = 1./(1+s*Lambda.^2);
        %--- RSS = Residual sum-of-squares
        if aow>0.9 % aow = 1 means that all of the data are equally weighted
            % very much faster: does not require any inverse DCT
            RSS = norm(DCTy(:).*(Gamma(:)-1))^2;
        else
            % take account of the weights to calculate RSS:
            yhat = idctn(Gamma.*DCTy);
            RSS = norm(sqrt(Wtot(IsFinite)).*(Raw(IsFinite)-yhat(IsFinite)))^2;
        end
        %---
        TrH = sum(Gamma(:));
        GCVscore = RSS/numOfFinites/(1-TrH/numOfElements)^2;
    end
    
end
function z = InitialGuess(y,I)
    %-- nearest neighbor interpolation (in case of missing values)
    if any(~I(:))
        [~,L] = bwdist(I);
        z = y;
        z(~I) = y(L(~I));
    else
        z = y;
    end
    %-- coarse fast smoothing using one-tenth of the DCT coefficients
    siz = size(z);
    z = dctn(z);
    for k = 1:ndims(z)
        z(ceil(siz(k)/10)+1:end,:) = 0;
        z = reshape(z,circshift(siz,[0 1-k]));
        z = shiftdim(z,1);
    end
    z = idctn(z);
end
function [y,w] = dctn(y,DIM,w)
    y = single(y);
    sizy = size(y);
    
    % Test DIM argument
    if ~exist('DIM','var'), DIM = []; end
    assert(~isempty(DIM) || ~isscalar(DIM),...
        'DIM must be a scalar or an empty array')
    assert(isempty(DIM) || DIM==round(DIM) && DIM>0,...
        'Dimension argument must be a positive integer scalar within indexing range.')
    
    if isempty(DIM), y = squeeze(y); end
    dimy = ndims(y);
    
    % Some modifications are required if Y is a vector
    if isvector(y)
        dimy = 1;
        if size(y,1)==1
            if DIM==1, w = []; return
            elseif DIM==2, DIM=1;
            end
            y = y.';
        elseif DIM==2, w = []; return
        end
    end
    
    % Weighting vectors
    if ~exist('w','var') || isempty(w)
        w = cell(1,dimy);
        for dim = 1:dimy
            if ~isempty(DIM) && dim~=DIM, continue, end
            n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
            w{dim} = exp(1i*(0:n-1)'*pi/2/n);
        end
    end
    
    % --- DCT algorithm ---
    if ~isreal(y)
        y = complex(dctn(real(y),DIM,w),dctn(imag(y),DIM,w));
    else
        for dim = 1:dimy
            if ~isempty(DIM) && dim~=DIM
                y = shiftdim(y,1);
                continue
            end
            siz = size(y);
            n = siz(1);
            y = y([1:2:n 2*floor(n/2):-2:2],:);
            y = reshape(y,n,[]);
            y = y*sqrt(2*n);
            y = ifft(y,[],1);
            y = bsxfun(@times,y,w{dim});
            y = real(y);
            y(1,:) = y(1,:)/sqrt(2);
            y = reshape(y,siz);
            y = shiftdim(y,1);
        end
    end
    y = reshape(y,sizy);
end

function [y,w] = idctn(y,DIM,w)
    y = single(y);
    sizy = size(y);
    
    % Test DIM argument
    if ~exist('DIM','var'), DIM = []; end
    assert(~isempty(DIM) || ~isscalar(DIM),...
        'DIM must be a scalar or an empty array')
    assert(isempty(DIM) || DIM==round(DIM) && DIM>0,...
        'Dimension argument must be a positive integer scalar within indexing range.')
    
    if isempty(DIM), y = squeeze(y); end
    dimy = ndims(y);
    
    % Some modifications are required if Y is a vector
    if isvector(y)
        dimy = 1;
        if size(y,1)==1
            if DIM==1, w = []; return
            elseif DIM==2, DIM=1;
            end
            y = y.';
        elseif DIM==2, w = []; return
        end
    end
    
    % Weighing vectors
    if ~exist('w','var') || isempty(w)
        w = cell(1,dimy);
        for dim = 1:dimy
            if ~isempty(DIM) && dim~=DIM, continue, end
            n = (dimy==1)*numel(y) + (dimy>1)*sizy(dim);
            w{dim} = exp(1i*(0:n-1)'*pi/2/n);
        end
    end
    
    % --- IDCT algorithm ---
    if ~isreal(y)
        y = complex(idctn(real(y),DIM,w),idctn(imag(y),DIM,w));
    else
        for dim = 1:dimy
            if ~isempty(DIM) && dim~=DIM
                y = shiftdim(y,1);
                continue
            end
            siz = size(y);
            n = siz(1);
            y = reshape(y,n,[]);
            y = bsxfun(@times,y,w{dim});
            y(1,:) = y(1,:)/sqrt(2);
            y = ifft(y,[],1);
            y = real(y*sqrt(2*n));
            I = (1:n)*0.5+0.5;
            I(2:2:end) = n-I(1:2:end-1)+1;
            y = y(I,:);
            y = reshape(y,siz);
            y = shiftdim(y,1);
        end
    end
    
    y = reshape(y,sizy);
end







function xy=findFFTPeak(I)
    siz1=size(I,1);
    siz2=size(I,2);
    I=sum(I(round(siz1/2-siz1*0.04):round(siz1/2+siz1*0.04),:),1);
    
    bwLabel=zeros(size(I));
    m=1; th=10:-0.1:0;
    while max(bwLabel)<2 && m<=numel(th)
        bw=I>(median(I)+std(I)*th(m));
        if mod(siz2,2)==0
            bw((siz2/2):end)=false;
        else
            bw(((siz2+1)/2):end)=false;
        end
        bwLabel=labelBW(bw);
        m=m+1;
    end
    
    if max(bwLabel)<2
        xy=[]; return
    else
        bw=bwLabel==max(bwLabel)-1;
        xy=find(I==max(I(bw)));
        xy=[round(siz1/2),xy(end)];
    end
end

function bwLabel=labelBW(bw)
    bwLabel=zeros(size(bw));
    if bw(1)
        m=1;
    else
        m=0;
    end
    bwLabel(1)=m;
    for k=2:numel(bw)
        if (bw(k-1) && bw(k))
            bwLabel(k)=m;
        elseif (~bw(k-1) && bw(k))
            m=m+1;
            bwLabel(k)=m;
        end
    end
end


