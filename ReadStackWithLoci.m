function varargout=ReadStackWithLoci(file,StackID,Range)
    % ReadStackWithLoci reads all kinds of file formats supported by
    % bio-formats.
    %
    % ====================================================================
    % In theory this function can read all kinds of file formats supported
    %   by bio-formats. Note that this function is not fully tested. Only
    %   Leica file format (*.lif) is tested.
    % ====================================================================
    % Also be informed that a LOCI_TOOLS.JAR file should be present in the
    %   same folder as this function. You can download a stable version
    %   from http://loci.wisc.edu/bio-formats/downloads.
    % ====================================================================
    %
    % Usage:
    %   1. Read single stack:   [I,Info]=ReadStack(file);
    %   2. Read specific stack: [I,Info]=ReadStack(file,StackID);
    %   3. Read specific frames in specific stack:
    %                           [I,Info]=ReadStack(file,StackID,[Start, Stop]);
    %   4. Get stack number :   StackNum=ReadStack(file,'GetStackNum');
    %   5. Get all meta data:   Metadata=ReadStack(file,'GetMetaData');
    % Output:
    %   Info=[SizeX, Pixel Physical Size X];
    %        [SizeY, Pixel Physical Size Y];
    %        [SizeT, Pixel Physical Size T];
    %        [SizeZ, Pixel Physical Size Z];
    %        [SizeC,                      ];
    
    
    %% Prepare outputs if returned in-btween.
    if nargout==1;      varargout{1}=[];
    elseif nargout==2;  varargout{1}=[]; varargout{2}=[];
    end
    
    if ((nargin==1) || (nargin==2) || (nargin==0));   Range=[];    end
    %% Check Bio-format library.
    % status = bfCheckJavaPath(1);
    
    jPath = javaclasspath('-all');
    isLociTools = cellfun(@(x) ~isempty(regexp(x, '.*bioformats_package.jar$', 'once')),jPath);
    status = any(isLociTools);
    if ~status
        % Assume the jar is in Matlab path or under the same folder as this file
        path = which('bioformats_package.jar');
        javaaddpath(path);  % Add loci_tools to dynamic Java class path
        loci.common.DebugTools.enableLogging('OFF'); % Options: DEBUG, OFF, INFO.
        status = true;
    end
    
    if ~status
        disp(['Missing Bio-Formats library. Either add bioformats_package.jar '...
            'to the static Java path or add it to the Matlab path.']);
        return
    end
    clear('isLociTools','jPath','status')
    
    %% Input.
    if nargin == 0
        [file, path] = uigetfile(TestedLociFormat, 'Choose a file to open');
        if isequal(path, 0) || isequal(file, 0), return; end
        file = fullfile(path, file);
    elseif exist(file, 'file') == 2
        
    elseif exist(file, 'dir') == 7
        [file, path] = uigetfile(TestedLociFormat, 'Choose a file to open',file);
        if isequal(path, 0) || isequal(file, 0), return; end
        file = fullfile(path, file);
    else
        [file, path] = uigetfile(TestedLociFormat, 'Choose a file to open');
        if isequal(path, 0) || isequal(file, 0), return; end
        file = fullfile(path, file);
    end
    
    if nargin<2; StackID=[]; end
    
    %% Open the file.
    % r=bfGetReader(file);
    
    r = loci.formats.ChannelFiller();
    r = loci.formats.ChannelSeparator(r);
    r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
    r.setId(file);
    
    
    %% Check if the image frame is too large.
    planeSize = loci.formats.FormatTools.getPlaneSize(r);
    if planeSize/(1024)^3 >= 2,
        error(['Image plane too large. Only 2GB of data can be extracted '...
            'at one time. You can workaround the problem by opening '...
            'the plane in tiles.']);
    end
    clear('planeSize')
    
    %% Get the image stack number.
    numSeries = r.getSeriesCount();
    
    %% Set pointer to image stack s, counting from 1,2,3,....
    if ischar(StackID)
        if strcmpi(StackID,'GetStackNum')
            if nargout>0; varargout=cat(1,{numSeries},cell(numel(2:nargout),1)); end
        elseif strcmpi(StackID,'GetMetaData')
            Metadata=char(r.getMetadata);
            Metadata=prepareMeta(Metadata);
            if nargout>0; varargout=cat(1,{Metadata},cell(numel(2:nargout),1)); end
        % elseif strcmpi(StackID,'GetXYZTC')
        end
        r.close();
        
        for k=2:nargout; varargout{k}=[]; end           %#ok<AGROW>
        return;
        
    elseif isempty(StackID)
        if numSeries>1; StackID=GetSeriesNumber(numSeries); else  StackID=1; end
        if isempty(StackID); r.close(); return; end
    elseif isscalar(StackID) && (~isempty(StackID)) && (StackID<1 || StackID>numSeries)
        r.close();
        disp(['    Wrong Stack ID. Should be between 1 ~ ' num2str(numSeries) '.']);
        return;
    end
    r.setSeries(StackID - 1);
    
    %% Retrieve Raw data.
    numImages = r.getImageCount();
    if isempty(Range)
        numStart=1; numStop=numImages;
    elseif numel(Range)==2
        if Range(1)>Range(2);Range=sort(Range,'ascend');end
        if (Range(1)>=1)&&(Range(1)<=numImages)
            numStart=Range(1);
        else
            numStart=1;
        end
        if (Range(2)>=1)&&(Range(2)<=numImages)
            numStop=Range(2);
        else
            numStop=numImages;
        end
    else
        disp('    Range error. Read full stack instead.')
        numStart=1; numStop=numImages;
    end
    
    pixelType = r.getPixelType();
    bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
    fp = loci.formats.FormatTools.isFloatingPoint(pixelType);
    % isSigned = loci.formats.FormatTools.isSigned(pixelType);
    bppMax = power(2, bpp * 8);
    little = r.isLittleEndian();
    X=r.getSizeX(); Y=r.getSizeY();
    
    plane = r.openBytes(numStart-1, 0, 0, X, Y);
    arr = loci.common.DataTools.makeDataArray(plane,bpp, fp, little);
    
    I=zeros(r.getSizeY,r.getSizeX,numStop-numStart+1,class(arr)); %#ok<ZEROLIKE>
    I(:,:,1) = reshape(arr, [X,Y])';
    for k = (numStart+1):numStop
        plane = r.openBytes(k-1, 0, 0, X, Y);
        arr = loci.common.DataTools.makeDataArray(plane,bpp, fp, little);
        I(:,:,k-numStart+1) = reshape(arr, [X,Y])';
    end
    
    %% Get the metadata.
    MetaData = r.getMetadataStore();
    a=zeros(numImages,1);
    try
        for kkk=1:numImages
            a(kkk)=MetaData.getPlaneDeltaT(StackID-1,kkk-1).floatValue;
        end
        a=median(diff(a));
    catch err %#ok<*NASGU>
        a=0;
    end
    
    try
        c=MetaData.getPixelsPhysicalSizeX(StackID-1).value;
        c=c.doubleValue;
    catch err;
        c=0;
    end
    try
        d=MetaData.getPixelsPhysicalSizeY(StackID-1).value;
        d=d.doubleValue;
    catch err;
        d=0;
    end
    try
        e=MetaData.getPixelsPhysicalSizeZ(StackID-1).value;
        e=e.doubleValue;
    catch err;
        e=0;
    end
    
    Info=([r.getSizeX,r.getSizeY,r.getSizeT,r.getSizeZ,r.getSizeC;c,d,a,e,0])';
    
    %% Close the file.
    r.close();
    
    %% Outputs.
    if nargout==1;      varargout{1}=I;
    elseif nargout==2;  varargout{1}=I;varargout{2}=Info;
    end
end






function fom=TestedLociFormat
    fom={'*.tif;*.tiff','Tiff file format (*.tif;*.tiff)';...
        '*.tf8','Big Tiff format (*.tf8)';...
        '*.lif','Leica LAS AF LIF (*.lif)';...
        '*.lsm','Zeiss LSM (Laser Scanning Microscope) (*.lsm)';...
        '*.vws','Till Photonics TillvisION (*.vws)'};
    
    All=char(fom{:,1});
    All=cat(2,All,repmat(';',[size(All,1),1]));
    All=regexprep((reshape(All',[1,numel(All)])),' ','');
    All={All,['All tested formats (',All,')']};
    fom=cat(1,All,fom);
    
    fom=cat(1,fom,{'*.*','All files (*.*)'});
end

function s=GetSeriesNumber(NumOfSeries)
    str=cell(NumOfSeries,1);
    for k=1:NumOfSeries; str{k}=num2str(k);end
    s=listdlg('Name','Series...','PromptString','Select a Series to read:',...
        'SelectionMode','single','ListString',str);
end

function meta=prepareMeta(strMeta)
    % Remove useless {}.
    if strMeta(1)=='{';   strMeta=strMeta(2:end);   end
    if strMeta(end)=='}'; strMeta=strMeta(1:end-1); end
    
    % Remove the space after cormma.
    strMeta = regexprep(strMeta,', ',',');
    
    % Find the cormma.
    bw=strfind(strMeta, ',');
    
    % Convert into cells.
    bw=[0,bw,numel(strMeta)+1];
    meta=cell((numel(bw)-1),1);
    for k=1:(numel(bw)-1)
        meta{k}=strMeta((bw(k)+1):bw(k+1)-1);
    end
    meta=sort(meta);
end