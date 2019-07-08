function varargout=SparkSites(SparkPos,SparkProperty,xyt_dim,varargin)
    %% Generate Spark Cluster.
    % SparkSitePos=[x1,x2,y1,y2,z1,z2,xc,yc,z_onset,SparkSiteID];
    % SparkSiteProperty=[FWHM,FWHM_R2,dAmp,Bgr,Tau,Sigma,R2,dAmp_dF/F0,DetctionMass,
    %                dF/F0_Mass,SparkSiteID];
    % Spark_SiteRelation=[SparkSiteNo, SparkID];
    %
    % Usage:
    %     [SparkSitePos,SparkSiteProperty,Spark_SiteRelation,SparkCluster]=SparkSites(...
    %         S.SparkPos,S.SparkProperty,S.xyt_dim,'LiveView',true);
    
    %% Input.
    p=inputParser;
    p.addParameter('LiveView',0,@(x)isscalar(x) && x>=0)
    p.addParameter('MinDistance',1.0,@(x)isscalar(x) && x>=0)       % In um now.
    parse(p, varargin{:});
    p=p.Results;
    MinSparkSiteDist=p.MinDistance/xyt_dim(1);
    clear('p','varargin')
    
    % PrintOutResults('I',S.xyt_dim,'',[],S.Gain,S.DetectionLimit,S.TauLim,S.FWHMLimit,...
    % S.CellMask,S.Creator,S.SparkSitePos,S.SparkSiteProperty,S.SparkFrequency);
    
    %% Spark information.
    % SparkPos=[x1,x2,y1,y2,z1,z2,xc,yc,z_onset,SparkID];
    % SparkProperty=[FWHM,FWHM_R2,dAmp,Bgr,Tau,Sigma,R2,dAmp_dF/F0,DetctionMass,
    %                dF/F0_Mass,SparkID];
    
    %% Clustering.
    SparkSpatioalPos=SparkPos(:,7:8);
    SparkPdist=pdist(SparkSpatioalPos,'euclidean');
    SparkLinage=linkage(SparkPdist,'single');
    
    % Plot the hierarchical tree.
    % figure;dendrogram(SparkLinage,0);
    % title('Hierarchical Tree of All Sparks');
    % xlabel('Spark ID');ylabel('Distance (in pixels)');

    %SparkLinage=inconsistent(SparkLinage);
    SparkCluster=cluster(SparkLinage,'cutoff',MinSparkSiteDist,'Criterion','distance','depth',round(size(SparkPos,1)/2));
    clear('MinSparkSiteDist','xyt_dim','SparkLinage','SparkPdist','SparkSpatioalPos')
    
    SparkCluster=cat(2,SparkPos,SparkProperty,SparkCluster);
    size2_SparkPos=size(SparkPos,2);
    % clear('SparkPos','SparkProperty')
    SparkCluster=sortrows(SparkCluster,size(SparkCluster,2));
    
    %% Extract information for every spark site.
    SparkSiteNum=max(SparkCluster(:,end));
    SparkPropertyDataWidth=size(SparkCluster,2);
    SparkSiteProperty=zeros(SparkSiteNum,SparkPropertyDataWidth);
    Spark_SiteRelation=struct('SparkSiteNo',[],'SparkID','');
    parfor k=1:SparkSiteNum
        bw=SparkCluster(:,end)==k; %#ok<*PFBNS>
        CurrSite=SparkCluster(bw,:);
        SparkSiteProperty(k,:)=median(CurrSite,1);
        Spark_SiteRelation(k).SparkSiteNo=k;
        Spark_SiteRelation(k).SparkID=num2str(CurrSite(:,end-1)');
    end
    Spark_SiteRelation=Spark_SiteRelation';
    clear('CurrSite','SparkSiteNum','bw','k');
    
    SparkSitePos=cat(2,SparkSiteProperty(:,1:(size2_SparkPos-1)),SparkSiteProperty(:,end));
    SparkSiteProperty=SparkSiteProperty(:,(size2_SparkPos+1):end);
    SparkSiteProperty(:,end-1)=[];
    
    %% Output
    varargout(1)={SparkSitePos};
    varargout(2)={SparkSiteProperty};
    varargout(3)={Spark_SiteRelation};
    varargout(4)={SparkCluster};
    
end
