function varargout=SparkDsptSimAllInOneV1_2(varargin)
% Spark simulation with decriptive model.
%
% This function use the key/value inputs.
%
% Usage:
% [spark,linescan]=SparkDsptSimAllInOneV1_2;
%   or
% [spark,linescan]=SparkDsptSimAllInOneV1_2(...,'Property',value,...);
%
%     'FWHMMax':        1.5                 % In micrometer.
%     'Mu':             20                  % In millisecond.
%     'Hill':           0.6                 % In millisecond, for Hill funciton.
%     'Tau':            7                   % In millisecond.
%     'TauFWHM':        1.5                 % In millisecond.
%     'SpR':            0.056               % In micrometer, spatial resolution.
%     'TpR':            0.1                 % In millisecond, temporal resolution.
%     'DecayLimit':     0.03
%     'LiveView':       false
%     'OutputSettings': false

%% Standard settings for internal control.
    % p.addParamValue('Mu',20,@(x)isscalar(x));                             % In millisecond.
    % p.addParamValue('Sigma',3,@(x)isscalar(x));                           % In millisecond.
    % p.addParamValue('Tau',60,@(x)isscalar(x));                            % In millisecond.
    % % Parameters for the FWHM function.
    % p.addParamValue('FWHMMax',1.5,@(x)isscalar(x));                       % In micrometer.
    % p.addParamValue('TauFWHM',15,@(x)isscalar(x));                        % In millisecond.
%% Input.
    p=inputParser;
        % Parameters for the HillExp function.
        p.addParameter('Mu',10,@(x)isscalar(x));                             % In millisecond.
        p.addParameter('Sigma',3,@(x)isscalar(x));                           % In millisecond.
        p.addParameter('Tau',40,@(x)isscalar(x));                            % In millisecond.

        % Parameters for the FWHM function.
        p.addParameter('FWHMMax',1.5,@(x)isscalar(x));                       % In micrometer.
        p.addParameter('TauFWHM',10,@(x)isscalar(x));                        % In millisecond.

        p.addParameter('SpR',0.056,@(x)isscalar(x));                         % In micrometer.
        p.addParameter('TpR',1,@(x)isscalar(x));                             % In millisecond.

        p.addParameter('DecayLimit',0.001,@(x)isscalar(x) && x<=0.01 && x>0);

        p.addParameter('LiveView',false,@(x)islogical(x));
        p.addParameter('OutputSettings',false,@(x)islogical(x));
        parse(p, varargin{:});
    p=p.Results;
    clear('varargin')

%% Input settings
    if p.OutputSettings
        fprintf('\tSingle spark properties\n');
        fprintf('\t%-30s%2.2f um\n','Max spark FWHM:',p.FWHMMax);
        fprintf('\t%-30s%2.2f ms\n','Exponential FWHM Tau:',p.TauFWHM);

        fprintf('\t%-30s%2.2f ms\n','Spark delay:',p.Mu);
        fprintf('\t%-30s%2.2f ms\n','Upstroke Sigma:',p.Sigma);
        fprintf('\t%-30s%2.2f ms\n','Decay Tau:',p.Tau);
        fprintf('\t%-30s%0.2f um x %0.2f um x %0.2f ms\n','Resolution (x,y,t):',p.SpR,p.SpR,p.TpR);
    end
%% 2D Gaussian surface on X and Y, and ExpExp on temporal dimension.

    % p(1), amplitude; p(2), mu; p(3), exponential tau; p(4), gaussian sigma.
    %GauConvExp=@(x,p) p(1)/2.*exp((p(4)^2+2*p(3)*(p(2)-x))/2/p(3)^2).*(1-(p(4)^2+p(3)*(p(2)-x))./
    %abs(p(4)^2+p(3)*(p(2)-x)).*erf(abs(p(4)^2+p(3)*(p(2)-x))/sqrt(2)/p(4)/p(3)));
    GauConvExpFWHM=@(FWHMMax,T,Mu,TauFWHM) FWHMMax*(1-exp((-(T-Mu)./TauFWHM))).*((1-exp((...
                -(T-Mu)./TauFWHM)))>0)+1e-9;
    GauConvExp=@(x,AmpMax,Mu,Tau,Sigma) AmpMax/2.*exp((Sigma^2+2*Tau*(Mu-x))/2/Tau^2).*(1-(Sigma^2+...
                Tau*(Mu-x))./abs(Sigma^2+Tau*(Mu-x)).*erf(abs(Sigma^2+Tau*(Mu-x))/sqrt(2)/Sigma/Tau));
    GauConvExpSpark=@(X,Y,T,FWHMMax,Mu,TauFWHM,AmpMax,Tau,Sigma) exp(-(X.^2+Y.^2)/2./(...
                GauConvExpFWHM(FWHMMax,T,Mu,TauFWHM)).^2).*GauConvExp(T,AmpMax,Mu+Sigma*2,Tau,Sigma);

%% Dimensions.
    tlen=-p.Tau*log(p.DecayLimit)+round(p.Mu);
    T=0:p.TpR:tlen;

    PositivePart=0:p.SpR:p.FWHMMax*3;   NegativePart=sort(PositivePart(2:end),'descend')*(-1);
    XYrange=[NegativePart PositivePart];

    [X,Y,T]=ndgrid(XYrange,XYrange,T);
    clear('tlen','PositivePart','NegativePart','XYrange')

%% Generate spark.
    spark=GauConvExpSpark(X,Y,T,p.FWHMMax,p.Mu,p.TauFWHM,1,p.Tau,p.Sigma);
    AmpTrace=GauConvExp(reshape(T(1,1,:),numel(T(1,1,:)),1),1,p.Mu,p.Tau,p.Sigma);
    FWHMTrace=GauConvExpFWHM(p.FWHMMax,reshape(T(1,1,:),numel(T(1,1,:)),1),p.Mu,p.TauFWHM);
    clear('X','Y','T','HillexpampSpark','Amp','FWHM')

%% Normalization.
    % max_pos_AmpTrace=find(AmpTrace==max(AmpTrace));max_pos_AmpTrace=max_pos_AmpTrace(end);
    % min_AmpTrace=min(AmpTrace(max_pos_AmpTrace:end));
    % AmpTrace=(AmpTrace-min_AmpTrace)/(max(AmpTrace(:))-min_AmpTrace);
    % AmpTrace=AmpTrace.*(AmpTrace>0);
    % spark=(spark-min_AmpTrace)/(max(spark(:))-min_AmpTrace);
    % spark=spark.*(spark>0);
    AmpTrace=AmpTrace/max(AmpTrace);
    spark=spark/max(spark(:));
%% Remove zeros in the matrix / Shrink the matrix.
    spark_positive_index=find(spark>0);
    [spark_x,spark_y,spark_t]=ind2sub(size(spark),spark_positive_index);
    spark=spark(min(spark_x):max(spark_x),min(spark_y):max(spark_y),min(spark_t):max(spark_t));

%% Generate linescan of the spark.
    if p.LiveView || nargout==2
        sparkCenter=floor((size(spark,1)+1)/2);
        linescan=spark(sparkCenter,:,:);    linescan=squeeze(linescan);
    end
%% Plotting.
    if p.LiveView
        T=1:numel(AmpTrace);    T=T*p.TpR;
        figure('Position',FigSize);
        subplot(3,1,1);
            plot(T,AmpTrace(:),'Marker','.');
            title('Central pixel transient');xlabel('Time (ms)');ylabel('Fluo (a.u.)');
            xlim([0 numel(AmpTrace)*p.TpR])
        subplot(3,1,2);
            plot(T,FWHMTrace(:),'Marker','.');
            title('Spark FWHM');xlabel('Time (ms)');ylabel('FWHM length (\mum)')
            xlim([0 numel(FWHMTrace)*p.TpR])
        subplot(3,1,3);imagesc(linescan);
            set(gca,'TickLength',[0 0],'XTickLabel',[],'YTickLabel',[]);
            title('Line Scan Mode of Spark');colormap(colorBlueToRed);
            %axis image;
    end

%% Output.
    if nargout==1
        varargout(1) = {spark};
    end
    if nargout==2
        varargout(1) = {spark};
        varargout(2) = {linescan};
    end
end


%% Figure size to screen size.
function Pos=FigSize
    Screen=get(0,'ScreenSize');
    FigHeight=round(Screen(4)*0.8);
    FigWidth=round(Screen(3)*0.2);
    Pos=[round((Screen(3)-FigWidth-100)*0.5) round(Screen(4)-FigHeight-100) FigWidth FigHeight];
end

function Cmap=colorBlueToRed
Cmap=[0 0 0.5625;0 0 0.58984375;0 0 0.6171875;0 0 0.64453125;0 0 0.671875;0 0 0.69921875;...
    0 0 0.7265625;0 0 0.75390625;0 0 0.78125;0 0 0.80859375;0 0 0.8359375;0 0 0.86328125;...
    0 0 0.890625;0 0 0.91796875;0 0 0.9453125;0 0 0.97265625;0 0 1;0 0.0625 1;0 0.125 1;...
    0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;...
    0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 0.9375;0.125 1 0.875;...
    0.1875 1 0.8125;0.25 1 0.75;0.3125 1 0.6875;0.375 1 0.625;0.4375 1 0.5625;0.5 1 0.5;...
    0.5625 1 0.4375;0.625 1 0.375;0.6875 1 0.3125;0.75 1 0.25;0.8125 1 0.1875;0.875 1 0.125;...
    0.9375 1 0.0625;1 1 0;1 0.933333337306976 0;1 0.866666674613953 0;1 0.800000011920929 0;...
    1 0.733333349227905 0;1 0.666666686534882 0;1 0.600000023841858 0;1 0.533333361148834 0;...
    1 0.466666668653488 0;1 0.400000005960464 0;1 0.333333343267441 0;1 0.266666680574417 0;...
    1 0.200000002980232 0;1 0.133333340287209 0;1 0.0666666701436043 0;1 0 0];
end
