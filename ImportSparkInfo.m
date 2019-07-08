function varargout=ImportSparkInfo(varargin)
    % Import analyzed spark information.
    % Usage: sparkInfo=ImportSparkInfo('Key','Value');
    %    Key                Value
    %    File               path_to_file.txt
    %    Info               true/false. Print out SparkInfo titles.
    %    TransformToCell    true/false. transform final format into cell format.
    %
    %   Output:
    %    [meta,metasimple,allinfo]=ImportSparkInfo(...);
    
    %% Inputs.
    Input=inputParser;
    Input.addParameter('File','',@(x) ischar(x) || isempty(x));
    Input.addParameter('Info',true,@(x) (islogical(x)) && (numel(x)==1));
    Input.addParameter('TransformToCell',false,@(x) (islogical(x)) && (numel(x)==1));
    parse(Input, varargin{:});
    Input=Input.Results;
    clear('varargin')
    
    if exist(Input.File,'file') ==2                   % If the string input is a file name.
        [Path,File,ext] = fileparts(Input.File);
        File=[File,ext];
        if isempty(Path); Path=pwd; end
        FullName = fullfile(Path, File);
    elseif exist(Input.File,'dir')==7                 % If the string input is a folder name.
        Path=Input.File;
        [File,Path] = uigetfile({'*.txt;*.asc'},'Please select a text/ASCII file:',Path);
        if File==0; return; end
        FullName = fullfile(Path, File);
    else
        [File,Path] = uigetfile({'*.txt;*.asc'},'Please select a text/ASCII file:');
        if File==0
            varargout(1)={''};
            return;
        else
            FullName=fullfile(Path, File);
        end
    end
    clear('ext','Path','File')
    
    %%
    fprintf('    Importing file:\n    %s',FullName);
    Info=dir(FullName);
    fprintf('\n\n    File size:\t%s KB (%0.0f Bytes).',thousandSep(Info.bytes),Info.bytes);
    
    fID=fopen(FullName,'r');
    InfoArray = fread(fID,Info.bytes,'*char');
    fclose(fID);
    InfoArray=InfoArray';
    
    lineNo=(1:numel(InfoArray));
    lineNo=lineNo(InfoArray==newline); %char(10)
    lineEnd=lineNo-1;
    lineStart=[1,lineNo(1:(numel(lineNo)-1))+1];
    
    lineStart=lineStart';
    lineEnd=lineEnd';
    clear('Info','fID','lineNo');
    
    
    %% Count cell number.
    cellNum=numel(strfind(InfoArray,'File: '));
    if cellNum<1
        varargout(1)={[]};
        varargout(2)={[]};
        varargout(3)={[]};
        fprintf('   Loading Error. No line started with "File: " found.\n')
        return;
    end
    fprintf('\n\n    Cell number in this file: %d.\n',cellNum);
    
    if Input.Info
        fprintf(['\n\n    Format for SparkInfo field or the 3rd output:\n    ' ...
            'cx(1) cy(2) ct(3) Mass(squm*ms,4) FWHM(um,5) AxisMax(6) ' ...
            'AxisMin(7) AmpdFF0(8)\n    Tau(9) UpstrokeSigma(ms,10) AmpInt(11) ' ...
            'BgrInt(12) DistToMem(um,13) Orient(degree,14)\n']);
        if nargout>=2
            fprintf(['\n    Second output format (Metainfo format):\n    CellName, RawStart, ' ...
                'RowStop, ItemNumber, CellArea, Frequency.\n\n']);
        end
    end
    
    %% Scane every line.
    files(cellNum,1).FileName=[];
    cellNum=0;
    k=0;
    while k<numel(lineStart)
        k=k+1;
        if (lineEnd(k)-lineStart(k)+1)<2; continue; end
        
        currStr=InfoArray(lineStart(k):lineEnd(k));
        initialStr=currStr(1:6);
        switch initialStr
            case 'File: '
                cellNum=cellNum+1;
                [files(cellNum).FileName,files(cellNum).FilePath]=getFileName(currStr(7:end));
            case 'Detect'
            case 'xyt_di'
                currStr=currStr(8:end);                             % currStr=currStr(8:numel(currStr)-6);
                bw=((currStr>=48)&(currStr<=57))...                % currStr(currStr=='x')='';
                    | (currStr==46) | (currStr==32) | (currStr==9);
                currStr=currStr(bw);                                % currStr(currStr=='u')='';
                % currStr(currStr=='m')='';
                xyt_dim=str2num(currStr); %#ok<ST2NM>
                files(cellNum).xyt_dim=xyt_dim(1:3);
                if numel(xyt_dim)>=4
                    files(cellNum).xyt_len=xyt_dim(4:end);
                end
            case 'Thresh'
                files(cellNum).Threshold=str2double(currStr(13:end));
            case 'RawGai'
                files(cellNum).Gain=str2num(currStr(8:end)); %#ok<ST2NM>
            case 'ScaleO'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                TempStr(TempStr==':')=',';
                files(cellNum).ScaleOfDetection=str2num(TempStr); %#ok<ST2NM>
            case 'MassLi'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                files(cellNum).MassLimit=str2num(TempStr); %#ok<ST2NM>
            case 'AmpLim'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                files(cellNum).AmpLimit=str2num(TempStr); %#ok<ST2NM>
            case 'Leadin'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                files(cellNum).LeadingTail=str2num(TempStr); %#ok<ST2NM>
            case 'TauLim'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                files(cellNum).TauLimit=str2num(TempStr); %#ok<ST2NM>
            case 'FWHMLi'
                TempStr=currStr((find(currStr=='[')+1):(end-1));
                files(cellNum).FWHMLimit=str2num(TempStr); %#ok<ST2NM>
            case 'MinSpa'
                files(cellNum).MinSparkDist=str2double(currStr(26:end));
            case 'PhaseC'
                files(cellNum).PhaseContrast=str2double(currStr(end))==1;
            case 'ImageT'
                files(cellNum).ImageTranslation=str2double(currStr(end))==1;
            case 'CellAr'
                files(cellNum).CellArea=str2double(currStr(16:end));
            case 'ItemNu'
                files(cellNum).ItemNumber=str2double(currStr(12:end));
                fprintf('    No.%3g\tspNum=%4g\t%s\n',cellNum,files(cellNum).ItemNumber,files(cellNum).FileName);
            case 'ItemFr'
                files(cellNum).ItemFrequency=str2double(currStr(25:end));
            case 'Progra'
                files(cellNum).ProgramVersion=currStr(17:end);
            case '    No'
                if strcmp(currStr(9:11),'Var')
                    continue;
                end
                C=zeros(files(cellNum).ItemNumber,14);
                for j=1:files(cellNum).ItemNumber
                    currStr=InfoArray((lineStart(k+j-1)+4):lineEnd(k+j-1));
                    currStr((currStr=='I')|(currStr=='(')|(currStr==')')|(currStr==',')|(currStr==char(9)))=' ';
                    currStr=regexprep(currStr,'No.','');
                    currNums=str2num(currStr);   %#ok<ST2NM>
                    C(j,1:size(currNums,2)-1)=currNums(2:end);
                end
                k=k+files(cellNum).ItemNumber;
                files(cellNum).SparkInfo=C;
        end
    end
    
    
    m=0;
    for k=1:cellNum
        if ~isempty(files(k).ItemNumber)
            m=m+files(k).ItemNumber;
        end
    end
    fprintf('\n    In summary, %d sparks imported.\n\n',m);
    
    
    %% Combine all sparks for further usages.
    if nargout>=2
        m=0;
        for k=1:cellNum
            if ~isempty(files(k).ItemNumber)
                m=m+files(k).ItemNumber;
            end
        end
    
        meta=cell(numel(files),8);
        C=zeros(m,14);
        n=0;
        for k=1:cellNum
            meta{k,1}=files(k).FileName;
            meta{k,2}=files(k).FilePath;
            if isempty(files(k).ItemNumber); continue; end
            m=n+1;                       meta{k,3}=m;
            n=m+files(k).ItemNumber-1;   meta{k,4}=n;
    
            meta{k,5}=files(k).ItemNumber;
            meta{k,6}=files(k).CellArea;
            meta{k,7}=files(k).ItemFrequency;
            meta{k,8}=files(k).SparkInfo;
            try
            C(m:n,:)=files(k).SparkInfo;
            catch
                disp(n)
            end
        end
        meta=cat(1,{'File','Path','StartRow','StopRow','ItemNumber','CellArea','ItemFrequency','SparkInfo'},meta);
    end
    
    %% Calculate Median values for each cell.
    files=sparkOrganize(files);
    
    %% Output
    if nargout==1
        varargout(1)={files};
    elseif nargout==2
        varargout(1)={files};
        varargout(2)={meta};
    elseif nargout==3
        varargout(1)={files};
        varargout(2)={meta};
        varargout(3)={C};
    end
end

function files=sparkOrganize(files)
    for k=1:numel(files)
        if ~isempty(files(k).SparkInfo)
            files(k).MassOfCellMedian=median(files(k).SparkInfo(:,4),1);
            files(k).FWHMOfCellMedian=median(files(k).SparkInfo(:,5),1);
            files(k).AxisMaxOfCellMedian=median(files(k).SparkInfo(:,6),1);
            files(k).AxisMinOfCellMedian=median(files(k).SparkInfo(:,7),1);
            files(k).AmpFFzeroOfCellMedian=median(files(k).SparkInfo(:,8),1);
            files(k).TauOfCellMedian=median(files(k).SparkInfo(:,9),1);
            files(k).UpstrokeOfCellSigmaMedian=median(files(k).SparkInfo(:,10),1);
            files(k).AmpIntOfCellMedian=median(files(k).SparkInfo(:,11),1);
            files(k).BgrIntOfCellMedian=median(files(k).SparkInfo(:,12),1);
            files(k).DistToOfCellMemMedian=median(files(k).SparkInfo(:,13),1);
        end
    end
end

function [FileName,FilePath]=getFileName(str)
    bw=[find(str=='/'),find(str=='\')];
    if numel(bw)>0
        FileName=str((max(bw)+1):end);
        FilePath=str(1:(max(bw)));
    else
        FileName='NotKnown';
        FilePath='NotKnown';
    end
end

function s=thousandSep(n)
    n=n/1024; % kb
    MB=n/1024;
    KB=round(mod(n,1024)*10)/10;
    if MB>=1
        s=[num2str(floor(MB)),',',num2str(KB)];
    else
        s=num2str(KB);
    end
end
