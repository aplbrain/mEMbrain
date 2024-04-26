function compute2dsegs(membranePath,mip,sections,outFolder,reduceMin,crop,ds)


%%% if probability is sparse, this would make the hmin computation run in
%%% small windows to avoid heavy computation on constant regions.
sparse_prob = 1;

% nag's P0 cereb.

% ds for 2dseg is an input

ds_membrane = 2;

DEBUG = 0;
sectionStopper = 0;

if (0)
    
  
    
    sectionTime=tic; % dce04
    reduce = 0.015;
    ds = 1;
    crop = 0;
    mip=1;
    membranePath = fullfile('.','../membranes',{'GT1'  'GT2_preGT1'  'GT2_preGT1Skel'});
    compute2dsegs( ...
        membranePath, ...
        mip,   0:0, '2dseg-GT1_2',reduce,crop,ds);
    sectionTime_elapsed=toc(sectionTime);
    

    sectionTime=tic; % dce04
    reduce = 0.007;
    ds = 1;
    crop = 0;
    mip=1;
    membranePath = fullfile('.','../membranes',{'GT1'  'GT2_preGT1'  'GT2_preGT1Skel'});
    compute2dsegs( ...
        membranePath, ...
        mip,   0:93, '2dseg-GT1_2',reduce,crop,ds);
    sectionTime_elapsed=toc(sectionTime);
    
     
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




outFolder = sprintf('%s_%g_ds%d_cr%d',outFolder,reduceMin,ds,crop);

overSegmentation = 0;
out = fullfile('./../2dseg',outFolder);

if overSegmentation
    out = [out '_oversegmentation'];
end

fmt = 'png' %%% read



patternTiles_read = 'sect_%06d_r%d_c%d' %%%
%patternTiles_write = 'sect_%06d_r%d_c%d' %%%

%%% sections
patternSection_read = 'Sect_%06d';
patternSection_write = 'Sect_%06d';


%%% mip zero tiling: column and row
colmin = 0;
colmax = 34-1;
rowmin = 0;
rowmax = 34-1;

%%% mip tiling;
mipcolmin = floor(colmin/2^mip);
mipcolmax = ceil(colmax/2^mip);
miprowmin = floor(rowmin/2^mip);
miprowmax = ceil(rowmax/2^mip);


%%%
%sections = minSections:maxSections;
%minSections
%maxSections


mipMemPath = fullfile(membranePath, sprintf('mip%d',mip));


rng(7)

colorsuint32 = uint32([0, randperm(2^24-1)]);
colorsuint = reshape(typecast(colorsuint32,'uint8'),4,[]);
colors = permute(colorsuint(1:3,:), [2 1]);
colors(1,:) = 0;


smooth='nosmooth'; %'smooth'



outMipPath = fullfile(out, sprintf('mip%d_%s_notiles',mip,smooth));
mkdir(outMipPath);

outMipTilePath =  fullfile(out, sprintf('mip%d',mip));
mkdir(outMipTilePath)

outNextMipTilePath =  fullfile(out, sprintf('mip%d',mip+1));
mkdir(outNextMipTilePath)

%mxID = zeros(1,numel(sections),'uint64');

tileSize = [1024 1024];

%Nworkers = 24;
%p2=parpool(Nworkers)
%opts = parforOptions(p2)
Nworkers = 12;
 
%parfor (section_index = 1:numel(sections), Nworkers)   %(section_index = 1:numel(sections), opts)    
for section_index = 1:numel(sections)      
    for ipause=1:30 
        pause(0.05)
        sprintf('.................... %d:%d ------.........',Nworkers,section_index); 
    end
    
    sectionID = sections(section_index)
    sectionPath = fullfile(mipMemPath, sprintf(patternSection_read,sectionID));
    
    if DEBUG && sectionID >= sectionStopper
        keyboard
    end
    
    %keyboard
    
    sectionProb_cell = cell(1,numel(sectionPath));
    for imembPath=1:numel(sectionPath)
        sectionProb_cell{imembPath} = readSection(sectionPath{imembPath},mipcolmin,mipcolmax,miprowmin,miprowmax,crop, ...
            tileSize,patternTiles_read,sectionID,fmt);
    end
    sectionProbMax = max(cat(3,sectionProb_cell{:}),[],3);
    sectionProbMin = min(cat(3,sectionProb_cell{:}),[],3);
    sectionProb = uint8(single(sectionProbMax)*0.7+single(sectionProbMin)*0.3);
     
    
    outSectionPath = fullfile(outMipPath, sprintf(patternSection_read,sectionID));
    mkdir(outSectionPath);
    
    
    %keyboard
    'done tiling'
    % sectionProb = 255-sectionProb;
    
    %%% if requires cleaning from the outside
    if (0)
        frame=zeros(size(sectionProb),class(sectionProb));
        frame([1 end],:) = 1;
        frame(:,[1 end]) = 1;
        [y,x]=ind2sub(size(frame),find(frame));
        mask = bwselect(sectionProb==255,x,y);
        [y,x]=ind2sub(size(frame),find(sectionProb==0 & imdilate(mask,ones(3,3))));
        mask2 = bwselect(sectionProb==0,x,y);
        sectionProb(mask2) = 255;
    end
    
    %
    
    
    if strcmp(smooth,'smooth')
        vloc_s = double(imfilter(sectionProb,fspecial('gaussian',[13 13],1)))/255;
        'done filtering'
    else
        vloc_s = double(sectionProb)/255;
    end
    
    
    if crop > 0
        'cropping'
        tic
        vloc_s_crop = vloc_s(crop+1:end-crop,crop+1:end-crop);
        toc
        'croped'
        
    else
        'copying'
        tic
        vloc_s_crop = vloc_s;
        toc
        'copied'
    end
    
    
    % s0 = vloc_s_crop(12314:12314+1024,12314:12314+1024);
    
    %%% use this to generate smaller segments that are more likely to be
    %%% compaible across Z when a large component is spilt into two
    %%% components. The idea is to break the large component into two
    %%% patches and generate over-segmentation.
    
    if overSegmentation
        
        %%% keyboard
        tic
        %f=imregionalmax(imhmax(bwdist(bwareaopen(vloc_s>0.3,100)),1));
        %impose=imimposemin(vloc_s,f);
        
        %vloc_s_d = vloc_s;
        %vloc_s_d(vloc_s>0.4)=inf;
        %%d=graydist(vloc_s_d,f,'quasi-euclidean');
        %%impose(d>0.3)=inf;
        %w=watershed(impose,8);
        %%w2=w;
        %%w2(d>0.6)=0;
        t=imhmax(bwdist(bwareaopen(vloc_s>0.3,100)),1);
        g=uint32(watershed(max(t(:))-t,8));
        g(vloc_s>0.25)=0;
        toc
        
        
        %%%figure; im(labeloverlay(vloc_s,w2)); axis(a1)
        
        
    else
        % keyboard
        
        
        %vloc_min_sup = zeros(size(vloc_s_crop),class(vloc_s_crop));
        % s0; tic; vloc_min_sup= imhmin(s0,reduceMin,8); toc
        %
        
        'computing imhmin...'
        tic;
        if sparse_prob
            'computing imhmin...'
            window = 4096;
            overlap = 256;
            vloc_min_sup=computeMinWindows(vloc_s_crop, window, overlap,reduceMin);
            
        else
            
            vloc_min_sup= ...
                imhmin(vloc_s_crop,reduceMin,8);
            
        end
        toc
        'done imhmin'
        
        %vloc_min_sup(vloc_s_crop>0.95) = inf; % was 0.7 before jan5
        
        %vloc_min_sup(vloc_s_crop>0.95) = 1; % was 0.7 before jan5
        tic;
        removeMask = vloc_s_crop>0.95;
        toc
        'mask computed'
        
        
        
        
        %'computing unique'
        %tic;
        %uq_vloc = unique(vloc_min_sup);
        %toc
        %'unique computed'
        
        %if numel(uq_vloc) > 1
        
        %%% this settings with other settings produced errors
        %window = 2048;
        %overlap = 256;
        
        %%% this settings with other settings produced less or no errors 
        window = 4096;
        overlap = 512;
            
        %tic; g = uint32(watershed(vloc_min_sup,8)); toc
        
        'computing watershed'
        t1=tic;
        g=computeWSWindows_external(vloc_min_sup, window, overlap, removeMask);
        toc(t1);
        'watershed computed'
        
         

        % keyboard
        
        if (0)
            g_db=computeWSWindows(vloc_min_sup(20000:24000,20000:24000), window, overlap, removeMask(20000:24000,20000:24000));  
            
            'computing watershed'
            t1=tic;
            window_db = 4096;
            overlap_db = 512;
            g=computeWSWindows(vloc_min_sup, window_db, overlap_db, removeMask);
            toc(t1);
            'watershed computed'
            
            imwrite(uint8(mod(g(20000:24000,20000:24000),256)),'./ws2_crop.png');     
        end
               
        %else
        %	g = zeros(size(vloc_min_sup),'uint32');
        
        %end
        
        %%% keyboard
        % g_debug = g;        
        t1=tic
        'removing high probabiliy border from watershed...'

        
        bigMxArea = imreconstruct(imerode(sectionProbMax>0.95*255 ,ones(51,51)),sectionProbMax>0.95*255);
        
        S = vloc_s_crop>0.71 | isinf(vloc_min_sup) | bigMxArea;
        g(S) = 0; % was 0.7 before jan5
        'done'
        toc(t1)
        
        ug = setdiff(g,0);
        if ~isempty(ug)
            tic; H=histcounts(g,[ug; inf]); toc;
            'done histcounts'
            g(ismember(g,ug(H<70)))=0; % was 50 before jan5
        end
    end
    
    
    
    %%% g(g>0) = g(g>0);
    
    if section_index == 1
        
        try
            %mxID(section_index) = max(g(:));
        catch
            keyboard
        end
        
    else
        
        %mxID(section_index) = mxID(section_index-1) + uint64(max(g(:)));
        
    end
    
    
    
    %%%t=permute(reshape(typecast(g(:),'uint16'),[4 size(g)]),[2 3 1]);
    
    % keyboard
    
    if crop == 0
        g_full = g;
    else
        % copying cropped part'
        tic
        g_full = zeros(size(vloc_s),class(g));
        g_full(crop+1:end-crop,crop+1:end-crop) = g;
        toc
        'cropped part copied'
    end
    
    tic
    if ds == 1
        g_ds = g_full;
    elseif ds == 2
        g_ds = g_full(1:2:end,1:2:end);
    elseif ds == 4
        g_ds = g_full(1:4:end,1:4:end);
    end
    toc
    
    
    %%%% saving using MAT files: good if need to avoid small files (SWP?)
    if (0)
        
        t1=tic
        fname=fullfile(outSectionPath, sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d_%s.mat',sectionID,reduceMin,mip,ds,smooth));
        save(fname,'g_ds','mxID','sections','-v7.3');
        toc(t1)
        
        %%% t1=tic; d=load(fname); toc(t1);
        
        
        %t_save = tic;
        %tic; g_dsuint8 = reshape(typecast(g_ds(:),'uint8'),[4 size(g_ds)]); toc;
        %%%tic; g_dsuint_RGB = permute(g_dsuint8(1:3,:,:),[2 3 1]); toc
        %tic; g_dsuint_RGB = cat(3,squeeze(g_dsuint8(1,:,:)),squeeze(g_dsuint8(2,:,:)),squeeze(g_dsuint8(3,:,:))); toc
        %t1=tic; writeSection_dim(g_dsuint_RGB,'./test/mip1/Sect_000060/',sectionID); toc(t1);
        %toc(t_save)
        
        %t1=tic; sectionData=readSection_dim('./test',mip,sectionID,[0 84],[0 104], [0 inf 0 inf],3); toc(t1);
        
        %reshape(g_dsuint8,[4 size(g_ds)])
        
        'saved mat file'
        
        'saving additional mip level...'
        tic
        if ds == 1
            
            if exists('var','mxID')
                t_save_mip2=tic;
                g_ds = g_ds(1:2:end,1:2:end);
                ds2 = 2;
                save(fullfile(outSectionPath, ...
                    sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d_%s.mat',sectionID,reduceMin,mip,ds2,smooth)),'g_ds','mxID','sections','-v7.3');
                toc(t_save_mip2);
            else
                t_save_mip2=tic;
                g_ds = g_ds(1:2:end,1:2:end);
                ds2 = 2;
                save(fullfile(outSectionPath, ...
                    sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d_%s.mat',sectionID,reduceMin,mip,ds2,smooth)),'g_ds','sections','-v7.3');
                toc(t_save_mip2);
            end
            
        end
        toc
        'saved'
        
    end
    
    %tic
    %imwrite(t(:,:,1:3),fullfile(outSectionPath, ...
    %    sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d_48bit.png',sectionID,reduceMin,mip,1)));
    %toc
    
    %tic
    %rgb = uint8(reshape(colors(g+1,:),[size(g) 3]));
    %toc
    
    
    if (0)
        try
            tic
            colorSection = uint8(reshape(colors(mod(g_full,size(colors,1))+1,:),[size(g_full) 3]));
            
            imwrite(colorSection(1:ds:end,1:ds:end,:),fullfile(outSectionPath,  sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d.png',sectionID,reduceMin,mip,ds)));
            toc
        catch
            keyboard
        end
        'saved 2d segs downsampled'
    end
    
    %tic
    %imwrite(rgb,fullfile(outSectionPath, ...
    %    sprintf('2dSeg_sect%d_reduce%g_mip%d_ds%d.png',sectionID,reduceMin,mip,1)));
    %toc
    
    
    % required if the composite membranes are needed for later computation
    t1=tic;
    imwrite(sectionProb(1:ds_membrane:end,1:ds_membrane:end),fullfile(outSectionPath,'./memb.png'));
    toc(t1)
    'saved membrane'
    
    
    
    
    %%% break into tiles for visualization in VAST as image layer.
    %%% this will map the IDs into unique random RGB colors, just for visualization
    %%% purposes if max ID is larger than 2^24 and for analysis if less
    %%% than 2^24 objects.
    
     
    sectionFolder = fullfile(outMipTilePath,sprintf('Sect_%06d',sectionID));
    mkdir(sectionFolder);
    
    
    % imwrite(uint8(mod(g_full(20000:24000,20000:24000),256)),'./g_full_crop.png');   
    
    t_tiling=tic;
    yi=0;
    for ystart=1:tileSize(1):size(g_full,1)
        yi
        xi = 0;
        for xstart=1:tileSize(2):size(g_full,2)
            xi;
            tile = zeros([tileSize 3],'uint8');
            yend = min(ystart+tileSize(1)-1,size(g_full,1));
            xend = min(xstart+tileSize(2)-1,size(g_full,2));
            
            
            
            objTile = uint8(reshape(colors(mod(g_full(ystart:yend,xstart:xend), ...
                size(colors,1))+1,:),[yend-ystart+1 ,xend-xstart+1 3]));
            
            %objTile = colorSection(ystart:yend,xstart:xend,:);
            
            
            tile(1:size(objTile,1),1:size(objTile,2),:)=objTile;
            
            if sum(tile(:)) == 0
                xi = xi + 1;
                continue
            end
            
            tilepath=fullfile(sectionFolder,sprintf('sect_%06d_r%d_c%d.png',sectionID,yi,xi));
            imwrite(tile,tilepath,'png');
            xi = xi + 1;
        end
        yi = yi + 1;
    end
    toc(t_tiling)
    
    %%%% tiling next mip level
    sectionFolder = fullfile(outNextMipTilePath,sprintf('Sect_%06d',sectionID));
    mkdir(sectionFolder);
    
    g_ds2=g_ds(1:2:end,1:2:end);
    
    g_ds2_color = uint8(reshape(colors(mod(g_ds2, ...
                size(colors,1))+1,:),[size(g_ds2) 3]));
    
    % imwrite(g_ds2_color(10000:12000,10000:12000,:),'./g_ds2_crop_color.png'); 
    % imwrite(uint8(mod(g_ds2(10000:12000,10000:12000),256)),'./g_ds2_crop.png'); 
    t_tiling=tic;
    yi=0;
    for ystart=1:tileSize(1):size(g_ds2,1)
        yi
        xi = 0;
        for xstart=1:tileSize(2):size(g_ds2,2)
            xi;
            tile = zeros([tileSize 3],'uint8');
            yend = min(ystart+tileSize(1)-1,size(g_ds2,1));
            xend = min(xstart+tileSize(2)-1,size(g_ds2,2));
            
            
            
            objTile = uint8(reshape(colors(mod(g_ds2(ystart:yend,xstart:xend), ...
                size(colors,1))+1,:),[yend-ystart+1 ,xend-xstart+1 3]));
            
            %objTile = colorSection(ystart:yend,xstart:xend,:);
            
            
            tile(1:size(objTile,1),1:size(objTile,2),:)=objTile;
            
            if sum(tile(:)) == 0
                xi = xi + 1;
                continue
            end
            
            tilepath=fullfile(sectionFolder,sprintf('sect_%06d_r%d_c%d.png',sectionID,yi,xi));
            imwrite(tile,tilepath,'png');
            xi = xi + 1;
        end
        yi = yi + 1;
    end
    toc(t_tiling)
    
    
    
    %t1_mip=tic; preparemipmap_path(out,mip+1,mip+1,sectionID,'nearest', 3,0,1024); toc(t1_mip);
    
    
    
end

function sectionProb = readSection(sectionPath,mipcolmin,mipcolmax,miprowmin,miprowmax,crop, ...
    tileSize,patternTiles_read,sectionID,fmt)

sectionProb = 255*ones((miprowmax+1)*tileSize(1),(mipcolmax+1)*tileSize(1),'uint8');

%%% read section
for col=mipcolmin:mipcolmax
    col;
    for row=miprowmin:miprowmax
        tilename = fullfile(sectionPath,sprintf([patternTiles_read '.' fmt],sectionID,row,col));
        
        if ~exist(tilename,'file')
            continue
        end
        
        if (row+1)*tileSize(1) < crop || (col+1)*tileSize(2) < crop ...
                || row*tileSize(1)+1 > size(sectionProb,1) - crop ...
                || col*tileSize(2)+1 > size(sectionProb,2) - crop
            continue
        end
        
        I=imread(tilename);
        
        sectionProb(row*tileSize(1)+1:(row+1)*tileSize(1), ...
            col*tileSize(2)+1:(col+1)*tileSize(2)) = I;
        
    end
end

function vol_hmin=computeMinWindows(vol, window, overlap,reduceMin)


vol_hmin = zeros(size(vol),class(vol));


for ystart=1:window:size(vol,1)
    
    
    for xstart=1:window:size(vol,2)
        
        xstartp = max(xstart-overlap,1);
        ystartp = max(ystart-overlap,1);
        xendp = min(xstart+window-1+overlap, size(vol,2));
        yendp = min(ystart+window-1+overlap, size(vol,1));
        
        patch = vol(ystartp:yendp,xstartp:xendp);
        
        mn = min(patch(:));
        if max(patch(:))==mn
            imhmin_patch = mn;
        else
            imhmin_patch=imhmin(patch,reduceMin,8);
        end
        imhmin_patch=imhmin_patch((ystart-ystartp)+1:end,(xstart-xstartp)+1:end);
        
        vol_hmin(ystart:ystart+size(imhmin_patch,1)-1,xstart:xstart+size(imhmin_patch,2)-1) ...
            = imhmin_patch;
        
    end
end


function vol_ws_cc=computeWSWindows(vol, window, overlap, removeMask)
% vloc_min_sup
%window = 2048;
%overlap = 256;

vol_ws1 = zeros(size(vol),'uint32');
mxID = uint32(0);

%checker = false(size(vol));

t1=tic;
yi=0;
for ystart=1:window:size(vol,1)
    yi = yi + 1;
    
    if mod(yi,2)==0
        x0 = 1;
    else
        x0 = 1+window;
    end
    
    ystart;
    for xstart=x0:2*window:size(vol,2)
        
        xstartp = max(xstart-overlap,1);
        ystartp = max(ystart-overlap,1);
        xendp = min(xstart+window-1+overlap, size(vol,2));
        yendp = min(ystart+window-1+overlap, size(vol,1));
        %xend = min(xstart+window-1, size(vol,2));
        %yend = min(ystart+window-1, size(vol,1));
        
        patch = vol(ystartp:yendp,xstartp:xendp);
        
        patchedMask = removeMask(ystartp:yendp,xstartp:xendp);
        
        
        if all(patchedMask,'all')
            %ws_patch = zeros(size(patch));
            continue
        else
            ws_patch=watershed(patch,8);
        end
        %ws_patch=ws_patch((ystart-ystartp)+1:end,(xstart-xstartp)+1:end);
        %box = true(size(patch,1)-(ystart-ystartp)-(yendp-yend), ...
        %    size(patch,2)-(xstart-xstartp)-(xendp-xend));
        %checker(ystart:ystart+size(box,1)-1,xstart:xstart+size(box,2)-1) = box;
        
        ws_patch_gl = uint32(ws_patch) + mxID;
        ws_patch_gl(ws_patch==0) = 0;
        
        %assert(max(ws_patch_gl(:)) == 0 || ...
        %  max(vol_ws1(:)) <  min(ws_patch_gl(ws_patch_gl>0)) );
        
        vol_ws1(ystartp:ystartp+size(ws_patch,1)-1,xstartp:xstartp+size(ws_patch,2)-1) ...
            = ws_patch_gl;
        mxID = mxID + max(uint32(ws_patch(:)));
        
        
        
    end
end
toc(t1)

vol_ws2 = zeros(size(vol),'uint32');

t1=tic;
yi=0;
for ystart=1:window:size(vol,1)
    yi = yi + 1;
    
    if mod(yi,2)==0
        x0 = 1+window;
        
    else
        x0 = 1;
    end
    
    ystart;
    for xstart=x0:2*window:size(vol,2)
        
        xstartp = max(xstart-overlap,1);
        ystartp = max(ystart-overlap,1);
        xendp = min(xstart+window-1+overlap, size(vol,2));
        yendp = min(ystart+window-1+overlap, size(vol,1));
        
        patch = vol(ystartp:yendp,xstartp:xendp);
        
        if isinf(min(patch(:)))
            %ws_patch = zeros(size(patch));
            continue
        else
            ws_patch=watershed(patch,8);
        end
        %ws_patch=ws_patch((ystart-ystartp)+1:end,(xstart-xstartp)+1:end);
        
        ws_patch_gl = uint32(ws_patch) + mxID;
        ws_patch_gl(ws_patch==0) = 0;
        
        %assert(max(ws_patch_gl(:)) == 0 || ...
        %  max(vol_ws1(:)) <  min(ws_patch_gl(ws_patch_gl>0)) );
        
        vol_ws2(ystartp:ystartp+size(ws_patch,1)-1,xstartp:xstartp+size(ws_patch,2)-1) ...
            = ws_patch_gl;
        mxID = mxID + max(uint32(ws_patch(:)));
        
        
        
    end
end
toc(t1)

t1=tic;
Q = vol_ws1>0 & vol_ws2;
pairs = unique([vol_ws1(Q) vol_ws2(Q)],'rows');
G=graph(pairs(:,1),pairs(:,2),1,max(vol_ws2(:)));
cc = conncomp(G);
cc = [0, cc];
toc(t1)

%keyboard
t1=tic;
vol_ws_cc=zeros(size(vol_ws1),'uint32');

yi=0;
for ystart=1:window:size(vol,1)
    yi = yi + 1;
    
    xi = mod(yi,2);
    
    ystart;
    for xstart=1:window:size(vol,2)
        xi = xi + 1;
        
        %xstartp = max(xstart-overlap,1);
        %ystartp = max(ystart-overlap,1);
        %xendp = min(xstart+window-1+overlap, size(vol,2));
        %yendp = min(ystart+window-1+overlap, size(vol,1));
        xend = min(xstart+window-1, size(vol,2));
        yend = min(ystart+window-1, size(vol,1));
        
        
        if mod(xi,2)~=0
            patch = vol_ws1(ystart:yend,xstart:xend);
            
        else
            patch = vol_ws2(ystart:yend,xstart:xend);
            
        end
        
        if max(patch(:)) > 0
            patch_cc = cc(patch+1);
        else
            continue
        end
        
        vol_ws_cc(ystart:yend,xstart:xend) = patch_cc;
        
        
    end
end
toc(t1)

%vol_ws=vol_ws2;
%vol_ws(checker) = vol_ws1(checker);
%vol_ws_cc = cc(vol_ws+1);







function [seg, allpairs]=matchWS(wsv,WeightThrehold,DemoMode)
% ws: path to h5
% Example:
% matchWS('/Users/yaronm/Dropbox (MIT)/MIT/yv_connectome/results/AC3/AC3-new-seeds-2D-WS/ws-new-supervoxels-AC3-TEST-NEW-SEEDS.h5',1.5);

% we define the bipartitle graph with binary edges (non-weighted) by
% weighing the similarity of segments and thresholdsing their normalized
% overlap.

allpairs = [];

seg = zeros(size(wsv),'uint32');

v1 = wsv(:,:,1);

for iz=1:size(wsv,3)-1
    
    v2 = wsv(:,:,iz+1);
    
    %%% VERIFY NO BUG -- after taking unique the objects IDs could potentially
    %%% be non-aligned between the images (this is OK for dense consecutive watershed as in the example)
    [obj1,~,IC1]  = unique(v1);
    [obj2,~,IC2]  = unique(v2);
    u1=uint16(reshape(IC1,size(v1)));
    u2=uint16(reshape(IC2,size(v2)));
    
    % defining the weight between each object in layer1 and layer2
    Doverlap1 = zeros(numel(obj1),numel(obj2));
    area1 = nan(numel(obj1),1);
    for i1=1:numel(obj1)
        
        area1(i1) = sum(u1(:) == i1);
        obj1lb_in_v2 = u2(u1 == i1);
        D_obj1_in_v2=histcounts(obj1lb_in_v2,[unique(obj1lb_in_v2); intmax('uint32')]);
        Doverlap1(i1, unique(obj1lb_in_v2)) = D_obj1_in_v2./area1(i1);
        
    end
    
    % do also the inverse and then symmetrize by
    % taking their sum
    Doverlap2 = zeros(numel(obj1),numel(obj2));
    area2 = nan(numel(obj2),1);
    for i2=1:numel(obj2)
        
        area2(i2) = sum(u2(:) == i2);
        obj2lb_in_v1 = u1(u2 == i2);
        D_obj2_in_v1=histcounts(obj2lb_in_v1,[unique(obj2lb_in_v1); intmax('uint32')]);
        Doverlap2(unique(obj2lb_in_v1),i2) = D_obj2_in_v1/area2(i2);
        
    end
    
    
    % maximization takes place here
    p6 = dmperm(Doverlap2+Doverlap1 > WeightThrehold); % 1.5
    
    %%%%%%%%
    pairs=unique([obj1(p6(p6~=0)), obj2(p6~=0)],'rows');
    pairs2=pairs(all(pairs,2),:);
    pairs2=unique(pairs2,'rows');
    
    allpairs = [allpairs; pairs2];
    %%%%%%
    
    
    
    % temporary; by-product of working with unique (compressed) values
    u2_t = p6(u2);
    I=u2_t==0;
    u2_t(I)=1;
    u2_t=obj1(u2_t);
    u2_t(I) = v2(I);
    
    if DemoMode
        
        M = [v1, u2_t];
        [y1,x1]=ind2sub(size(v1),find(imdilate(v1,ones(3))~=imerode(v1,ones(3))));
        [y2,x2]=ind2sub(size(v2),find(imdilate(v2,ones(3))~=imerode(v2,ones(3))));
        L=label2rgb(M,'hsv','k','shuffle');
        figure;
        subplot(121);  im(L(:,1:end/2,:)); colorbar off; hold on; plot(x1,y1,'.w');
        subplot(122);  im(L(:,end/2+1:end,:)); colorbar off; hold on; plot(x2,y2,'.w');
        set(gca,'fontsize',18);
        set(gca,'box','off')
        title('Lower Block Vs. Upper Block')
        
    end
    
    seg(:,:,iz) = v1;
    v1 = u2_t;
    
end
seg(:,:,end) = u2_t;

if DemoMode
    for iz=1:size(seg,3)
        imwrite(uint16(seg(:,:,iz)),fullfile(labelPath,strcat(num2str(iz,'%.4d'),'.tif')),'Tiff','compression','lzw');
    end
    
    scrollSegmentation(labelPath,'/Users/yaronm/Dropbox (MIT)/MIT/yv_connectome/ForestEdgeDetection/edges - Unix/ac3/ac3test');
    lb='/Users/yaronm/Dropbox (MIT)/MIT/yv_connectome/results/AC3/AC3_fully_conv_49x49_w0_dist_4_4M_3D_32f/test_3D/pipeline/stacked-labels-49x49-w0-dist-4-4M-ac3-3D-only-0.tif'
    
    lbv=readvolume(lb);
    [vi,vis,vim]=VI3D(seg,lbv); %#ok<*ASGLU>
    
    [vi,vis,vim]=VI3D(wsv,lbv);
    
    [vi,vis,vim]=VI3D(h5f,lbv);
    [vi,vis,vim]=VI3D(h5f,ws);
    supvox=readvolume(h5f);
    
    k=1;
    M = [supvox(:,:,k), supvox(:,:,k+1)];
    [y1,x1]=ind2sub(size(v1),find(imdilate(M(:,1:end/2,:),ones(3))~=imerode(M(:,1:end/2,:),ones(3))));
    [y2,x2]=ind2sub(size(v2),find(imdilate(M(:,end/2+1:end,:),ones(3))~=imerode(M(:,end/2+1:end,:),ones(3))));
    L=label2rgb(M,'hsv','k','shuffle');
    figure;
    subplot(121);  im(L(:,1:end/2,:)); colorbar off; hold on; plot(x1,y1,'.w');
    subplot(122);  im(L(:,end/2+1:end,:)); colorbar off; hold on; plot(x2,y2,'.w');
    set(gca,'fontsize',18);
    set(gca,'box','off')
    title('Lower Block Vs. Upper Block')
end






