function varargout=PhaseShiftAutoAlign(I,nthDim,xyshiftFit)
    %% PhaseShiftAutoAlign calculates and corrects phase shift generated in bi-directional scanning.
    % [I_corrected,xyshift]=ImAutoRegister(I,nthDim);
    
    %%
    I=single(I);
    
    if nargin==1
        nthDim=1;
    end
    if (nthDim~=1) && (nthDim~=2)
        disp('    dimension input is wrong.');
        return;
    end
    
    %% Calculate shift size based on the edge effects.
    if nargin<=2 || isempty(xyshiftFit)
        % Seperate into small average images.
        Image=mean(I,3);
        if nthDim==1
            Image1=Image(1:2:end,:);
            Image2=Image(2:2:end,:);
            if size(Image1,1)>size(Image2,1)
                Image1=Image1(1:(size(Image1,1)-1),:);
            end
            % Image1=diff(Image1,1,2);
            % Image2=diff(Image2,1,2);
        elseif nthDim==2
            Image1=Image(:,1:2:end);
            Image2=Image(:,2:2:end);
            if size(Image1,2)>size(Image2,2)
                Image1=Image1(:,1:(size(Image1,2)-1));
            end
            % Image1=diff(Image1,1,1);
            % Image2=diff(Image2,1,1);
        end
        
        % Translation calculation.
        ImDiff=@(I1,I2,p) fun_imDiff(I1,I2,p);
        options=optimset('MaxIter',100000000,'Display','off','TolX',1e-4);
        xyshift=[size(Image,1)*0.15,size(Image,2)*0.15];
        xyshiftFit=fminsearchbnd(@(p)ImDiff(Image1,Image2,p),[0,0],-xyshift,xyshift,options);
        clear('Image','Image1','Image2','xyshift');
    end
    %% Shift entire image;
    if nthDim==1
        xyshiftFit=xyshiftFit(2);
        for k=1:size(I,3)
            I(2:2:end,:,k)=ShiftImage(I(2:2:end,:,k),0,xyshiftFit);
        end
    elseif nthDim==2
        xyshiftFit=xyshiftFit(1);
        for k=1:size(I,3)
            I(:,2:2:end,k)=ShiftImage(I(:,2:2:end,k),0,xyshiftFit);
        end
    end
    
    %% Output.
    if nargout==1
        varargout{1}=I;
    elseif nargout==2
        varargout{1}=I;
        varargout{2}=xyshiftFit;
    end
end

function [R,I2]=fun_imDiff(I1,I2,p)
    I2=ShiftImage(I2,p(1),p(2));
    R=I1-I2;
    R=mean(R(:).^2);
end

function[Image] = ShiftImage(Image, deltax, deltay)
    % 01/25/99  Implemented by Edward Brian Welch, edwardbrianwelch@yahoo.com
    
    %% Pad the image.
    [oxdim,oydim,ozdim]=size(Image);
    deltax=deltax/oxdim; deltay=deltay/oydim;
    
    ImClass=class(Image);
    
    if mod(oxdim,2)==1 && mod(oydim,2)==1,
        tmp=zeros(oxdim+1,oydim+1,ozdim,ImClass);
        tmp(1:oxdim,1:oydim,:)=Image;
        Image=tmp;
    elseif mod(oxdim,2)==1 && mod(oydim,2)==0,
        tmp=zeros(oxdim+1,oydim,ozdim,ImClass);
        tmp(1:oxdim,1:oydim,:)=Image;
        Image=tmp;
    elseif mod(oxdim,2)==0 && mod(oydim,2)==1,
        tmp=zeros(oxdim+1,oydim,ozdim,ImClass);
        tmp(1:oxdim,1:oydim,:)=Image;
        Image=tmp;
    else
    end
    clear('tmp')
    
    % Put deltas into the range [0,2].
    while deltax<0,
        deltax = deltax + 2;
    end
    
    while deltax>2,
        deltax = deltax - 2;
    end
    
    while deltay<0,
        deltay = deltay + 2;
    end
    
    while deltay>2,
        deltay = deltay - 2;
    end
    
    if deltax~=0;Image=subshiftx(Image,deltax);end
    if deltay~=0;Image=subshifty(Image,deltay);end
    Image=real(Image);
    
    % Return an Image of the original size (in case it was odd)
    Image=Image(1:oxdim,1:oydim,:);
end

function Image=subshiftx(Image,deltax)
    xdim=size(Image,1);
    ydim=size(Image,2);
    if deltax==0; return; end
    % Calculate image's center coordinates
    xno = (xdim-1)/2;
    
    % Forward FFT image rows
    Image = fft(Image, xdim, 1);
    
    % Initialize constant part of the exponent expression.
    cons1 = (-2.0*pi/(xdim)) * deltax * xdim;
    
    % Calculate k values (Nyquist is at x=xno)
    xno_ceil=ceil(xno);
    k_array = zeros(xdim,1);
    k_array(1:xno_ceil)=0:xno_ceil-1;
    k_array((xno_ceil+1):xdim)=(xno_ceil-xdim):(-1);
    
    % Rotate the complex numbers by those angles
    angle_array = cons1*k_array;
    sin_ang = repmat(sin(angle_array),[1,ydim]);
    cos_ang = repmat(cos(angle_array),[1,ydim]);
    newr = real(Image).*cos_ang - imag(Image).*sin_ang;
    newi = real(Image).*sin_ang + imag(Image).*cos_ang;
    Image = newr + newi*1i;
    Image = ifft(Image, xdim, 1);
end

function Image=subshifty(Image,deltay)
    xdim=size(Image,1);
    ydim=size(Image,2);
    if deltay==0; return; end
    yno = (ydim-1)/2;
    
    Image =  fft(Image, ydim, 2);
    
    % Initialize constant part of the exponent expression.
    cons1 = (-2.0*pi/(ydim)) * deltay * ydim;
    
    yno_ceil=ceil(yno);
    k_array = zeros(ydim,1);
    k_array(1:yno_ceil)=0:yno_ceil-1;
    k_array((yno_ceil+1):ydim)=(yno_ceil-ydim):(-1);
    
    angle_array = cons1*k_array;
    sin_ang = repmat(sin(angle_array)',[xdim,1]);
    cos_ang = repmat(cos(angle_array)',[xdim,1]);
    newr = real(Image).*cos_ang - imag(Image).*sin_ang;
    newi = real(Image).*sin_ang + imag(Image).*cos_ang;
    Image = newr + newi*1i;
    Image = ifft(Image, ydim, 2);
end