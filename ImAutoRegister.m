function varargout=ImAutoRegister(I)
    %% All the frame in the I will be aligned to the first frame. 
    %   Intensity-based registration using imgegister in Matlab IPT.
    %   Syntax: I=ImAutoRegister(I);
    %% Registration.

    % fprintf('  0%%');
    [optimizer, metric] = imregconfig('monomodal'); % Intensity-based registration
    I_mean=mean(I,3);
    parfor k=1:size(I,3)
        % I(:,:,k)=imregister(I(:,:,k),I_mean,'translation',optimizer,metric);
        I(:,:,k)=imregister(I(:,:,k),I_mean,'affine',optimizer,metric);
        % fprintf('\b\b\b\b%3.0f%%',floor(k/TotalNum*100));
    end
    
    if nargout==1
        varargout{1}=I;
    end
end