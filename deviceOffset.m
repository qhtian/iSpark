function [Offset,Gain]=deviceOffset(I,Msize)
    %% Estimate device/detector gain value.
    %
    % Symtax: [Offset,Gain]=deviceGain(Raw,Msize);
    %   I:            Image stack need to be analyzed.
    %   MSize:        Optional. The grid size. Default is 8.
    
    %% New and faster processing.
    if nargin==1
        Msize=8;
    end
    % addpath('./mexFunctions');
    
    I=single(I);

    [x,y]=mxDeviceOffsetRearrange(I,Msize);
    
    x_sort=sort(x(:),'ascend');
    Offset=mode(x_sort(1:round(numel(x_sort)*0.2)));
    Offset=double(Offset);
    
    % Gain = x\y;
    params=wlsFit(x,y);
    Gain=params(1);
end

function params=wlsFit(x,y)
    params=zeros(2,1);
    IterativesMax=100;
    N=numel(x);
    Tol=1e-2;
    a0=10;    b0=-10;
    
    w=ones(N,1);
    for k=1:IterativesMax
        w2=w.^2;
        xy=w2.*x.*y;
        sw2 = sum(w2);
        sw2x=sum(w2.*x);
        sw2y=sum(w2.*y);
        params(1)=(sw2*sum(xy)-sw2x*sw2y);
        xy=w2.*x.^2;
        aux=sw2*sum(xy)-sw2x*sw2x;
        if abs(aux)<eps
            params(1)=a0;
        else
            params(1)=params(1)/aux;
        end
        
        params(2)=sw2y-params(1)*sw2x;
        if abs(aux)<eps
            params(2)=b0;
        else
            params(2)=params(2)/sw2;
        end
        r=y-x*params(1)-params(2);
        e=mean(abs(r));
        if (((abs(a0-params(1))<Tol)&&(abs(b0-params(2))<Tol)) || (e<Tol))
            break;
        end
        a0=params(1);
        b0=params(2);
        w=weightFun(r);
    end
end

function w=weightFun(x) % "huber" weighting function.
    x=x/InterQuarter(x);
    p=0.75;
    w=ones(size(x));
    bw=abs(x)>=p;
    if sum(bw)>0
        w(bw)=p./abs(x(bw)+eps);
    end
end

function range=InterQuarter(x)
    x=sort(x,'ascend');
    N=numel(x);
    m=floor(N/4);
    if m==0; m=1; end
    range=x(N-m)-x(m)+eps;
end