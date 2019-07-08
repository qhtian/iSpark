function [Iout,alpha,delta,sigma]=PUREDenoiseJava(I,varargin)
    %% This function needs PUREDenoise.jar file from BIG Imaging Center.
    %
    % Make sure that matlab can find both ij.jar and PureDenoise_jar.
    %
    % To add a jar to matlab:
    %	cd(prefdir)
    %	edit javaclasspath.txt
    %	javaclasspath % Display all the java class paths.
    %
    %	-- OR --
    %	edit classpath.txt % and put it in this file. It support $matlabroot. Use # to comment something.
    %
    %
    % Use Class in Java JAR File on Static Class Path
    %	This example shows how to call a class method in a Java?? Archive (JAR) file. The example uses a JAR file named mylibrary.jar in the folder C:\Documents\MATLAB\, containing a method, package.class.mymethod(params). Substitute your own JAR file name and path, and read the documentation to call your own method.
    %
    %	This example puts the JAR file on the static Java class path, making the classes always available in MATLAB??. If you only want to use the classes in the current MATLAB session, add the JAR file to the dynamic class path using the javaaddpath function.
    %
    %	Copy the mylibrary.jar file to the C:\Documents\MATLAB\ folder.
    %
    %	Add the JAR file to the static class path. Open the javaclasspath.txt file.
    %
    %		cd(prefdir)
    %		edit javaclasspath.txt
    %		Add the following text on a new line in the file.
    %			C:\Documents\MATLAB\mylibrary.jar
    %		Save and close the file. From now on, the JAR file is available to MATLAB.
    %		Restart MATLAB.
    %
    %	Call the method.
    %	package.class.mymethod(params)
    
    
    In=inputParser;
    In.addParameter('CS',6, @(x)(numel(x)==1 && x>0));      % Cycle of spins during desnoising.
    parse(In, varargin{:});
    In=In.Results;
    
    CS = In.CS;
    I=single(I);
    
    [Iout,alpha,delta,sigma]=PUREDenoiseJavaSingle(I,CS);
end

function [I,alpha,delta,sigma]=PUREDenoiseJavaSingle(I,CS)
    NBFRAME=1;
    Nmin=16;
    LOG=false;
    
    Isiz=[size(I,1),size(I,2),size(I,3)];
    if Isiz(1)<16 || Isiz(2)<16 || Isiz(3)<16
        return;
    end
    
    
    % denoise on x/y space.
    nxe = ceil(Isiz(1) / Nmin) * Nmin;
    nye = ceil(Isiz(2) / Nmin) * Nmin;
    Ext = [floor((nxe-Isiz(1))/2),floor((nye-Isiz(2))/2)];
    

    alpha=zeros(Isiz(3),1);
    delta=zeros(Isiz(3),1);
    sigma=zeros(Isiz(3),1);
    parfor k=1:Isiz(3)
        [I(:,:,k),alpha(k),delta(k),sigma(k)]=PURE_Denoise_Java_SingleFrame_Modul(I(:,:,k),Isiz,nxe,nye,CS,NBFRAME,LOG,Ext,0,0,0);
    end
    alpha=median(alpha);
    delta=median(delta);
    sigma=median(sigma);
    
    
    % denoise on x/t space
    I=permute(I,[1,3,2]);
    Isiz=[size(I,1),size(I,2),size(I,3)];
    nxe = ceil(Isiz(1) / Nmin) * Nmin;
    nye = ceil(Isiz(2) / Nmin) * Nmin;
    Ext = [floor((nxe-Isiz(1))/2),floor((nye-Isiz(2))/2)];
    parfor k=1:Isiz(3)
        [I(:,:,k),~,~,~]=PURE_Denoise_Java_SingleFrame_Modul(I(:,:,k),Isiz,nxe,nye,CS,NBFRAME,LOG,Ext,alpha,delta,sigma);
    end
    I=permute(I,[1,3,2]);
    
    % denoise on y/t space
    I=permute(I,[2,3,1]);
    Isiz=[size(I,1),size(I,2),size(I,3)];
    nxe = ceil(Isiz(1) / Nmin) * Nmin;
    nye = ceil(Isiz(2) / Nmin) * Nmin;
    Ext = [floor((nxe-Isiz(1))/2),floor((nye-Isiz(2))/2)];
    parfor k=1:Isiz(3)
        [I(:,:,k),~,~,~]=PURE_Denoise_Java_SingleFrame_Modul(I(:,:,k),Isiz,nxe,nye,CS,NBFRAME,LOG,Ext,alpha,delta,sigma);
    end
    I=permute(I,[3,1,2]);    
end

function [output,alpha,delta,sigma]=PURE_Denoise_Java_SingleFrame_Modul(I,Isiz,nxe,nye,CS,NBFRAME,LOG,Ext,alpha,delta,sigma)
    PUREdenoisingI=imageware.Builder.create(I,4);
    if (nxe ~= Isiz(1)) || (nye ~= Isiz(2))
        PUREdenoisingI = denoise.Operations.symextend2D(PUREdenoisingI, nxe, nye, Ext);
    end
    
    denoising = denoise.Denoising(PUREdenoisingI, alpha, delta, sigma, LOG, CS, NBFRAME);
    clear('PUREdenoisingI');
    % if alpha==0
    denoising.estimateNoiseParameters();
    % end
    denoising.perform();
    
    output = denoising.getOutput();
    output = output.getSliceDouble(0);
    output = reshape(output,[nxe,nye,1]);
    output = output((Ext(1)+1):(Ext(1)+Isiz(1)),(Ext(2)+1):(Ext(2)+Isiz(2)));
    output = single(output);
    
    alpha=denoising.getAlpha();
    delta=denoising.getDelta();
    sigma=denoising.getSigma();
    clear('denoising')
end










% function [I,alpha,delta,sigma]=PUREDenoiseJava3D(I,CS,NBFRAME)
%     Nmin=16;
%     LOG=false;
%     
%     alpha=zeros(size(I,3),1);
%     delta=zeros(size(I,3),1);
%     sigma=zeros(size(I,3),1);
%     
%     Isiz=[size(I,1),size(I,2),size(I,3)];
%     if Isiz(1)<16 || Isiz(2)<16
%         return;
%     end
% 
%     nxe = ceil(Isiz(1) / Nmin) * Nmin;
%     nye = ceil(Isiz(2) / Nmin) * Nmin;
%     Ext = [floor((nxe-Isiz(1))/2),floor((nye-Isiz(2))/2)];
% 
%     PUREdenoisingI=imageware.Builder.create(I,4);
%     if (nxe ~= Isiz(1)) || (nye ~= Isiz(2))
%         PUREdenoisingI = denoise.Operations.symextend2D(PUREdenoisingI, nxe, nye, Ext);
%     end
%     
%     denoising = denoise.Denoising(PUREdenoisingI, alpha, delta, sigma, LOG, CS, NBFRAME);
%     clear('PUREdenoisingI');
%     
%     denoising.estimateNoiseParameters();
% 
%     denoising.perform();
%     output = denoising.getOutput();
% 
%     for k=1:Isiz(3)
%         currSlice = output.getSliceDouble(k-1);
%         currSlice = reshape(currSlice,[nxe,nye,1]);
%         I(:,:,k) = currSlice((Ext(1)+1):(Ext(1)+Isiz(1)),(Ext(2)+1):(Ext(2)+Isiz(2)));
%     end
%     
%     alpha=median(denoising.getAlpha());
%     delta=median(denoising.getDelta());
%     sigma=median(denoising.getSigma());
%     clear('denoising','output');
% end