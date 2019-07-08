function SparkLabelRaw=CaEventSegmentWatershed(MaximaMask,xyt_dim)
    conn=26;
    sigma=0.42/xyt_dim(1);
    hsize=ceil(sigma*5);
    if mod(hsize,2)==0
        hsize=hsize+1;
    end
    h=fspecial('Gauss',hsize,sigma);
   
    I=-convn(single(MaximaMask),h,'same');
    SparkLabelRaw=watershed(I,conn);
    SparkLabelRaw(I>=-max(h(:)))=0;
end
