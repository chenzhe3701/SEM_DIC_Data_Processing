% estimate the grip direction based on xyTrans and xyStage info.
%
% chenzhe, 2017-05-22
%
% chenzhe, 2017-09-20 added note
% (1) This is truly a historic code for  20170430 test.
% (2) It looks at the hdr file for the four corner fovs (could look at all
% fovs).  Analyze the working distance, and stage position change.
% (3) The plot is a stage position change, indicating the axis of sample
% elongation -- this sample elongated along a direction 2.74 degrees ccw
% wrt the horizontal direction.
% (4) You could also look at WD, but looks like this sample was OK.

pathHdr = uigetdir('D:\Marissa_test_20170430\distortion_correction_setup_files\all hdr','choose the hdr folder');
pathDic = uigetdir('D:\Marissa_test_20170430\2]_20170430_ts5Al_02 tensile test_non_incremental_DIC','choose the dic folder');

if ~strcmpi(pathHdr(end),'\')
    pathHdr = [pathHdr,'\'];
end
if ~strcmpi(pathDic(end),'\')
    pathDic = [pathDic,'\'];
end

f1 = '20170430_ts5Al_02_e';
nE = 6;
nR = 3;
nC = 13;
micronPerPixel = 500/6144;
for ie = 0:6
    for ir = [0:3]
        for ic = [0,13]
            iE = ie + 1;
            iR = ir + 1;
            iC = ic + 1;
            
            fName = [f1,num2str(ie),'_r',num2str(ir),'c',num2str(ic)];
            fNameHdr = [fName,'.hdr'];
            fNameDic = [fName,'.mat'];
            
            WD(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'WD');
            StageX(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'StageX');
            StageY(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'StageY');
            
            WD(iR,iC,iE) = WD(iR,iC,iE) * 1e6;  % convert to micron
            StageX(iR,iC,iE) = StageX(iR,iC,iE) * 1e6; % make it to micron
            StageY(iR,iC,iE) = -StageY(iR,iC,iE) * 1e6;    % make it the same as image coordinate
            
            load([pathDic,fNameDic],'u','v');
            [nR,nC] = size(u);
            tR = round(nR/2);
            tC = round(nC/2);
            ImageU(iR,iC,iE) = u(tR,tC)*micronPerPixel;
            ImageV(iR,iC,iE) = v(tR,tC)*micronPerPixel;
            
            dW(iR,iC,iE) = WD(iR,iC,iE) - WD(iR,iC,1);
            dX(iR,iC,iE) = StageX(iR,iC,iE) - StageX(iR,iC,1) + ImageU(iR,iC,iE) - ImageU(iR,iC,1);
            dY(iR,iC,iE) = StageY(iR,iC,iE) - StageY(iR,iC,1) + ImageV(iR,iC,iE) - ImageV(iR,iC,1);
        end
    end
end
%%
figure; hold on;set(gca,'ydir','reverse');axis equal;
colors = jet(7);
data = [];
for ie = 0:6
    for ir = [0:3]
        for ic = [0,13]
            iE = ie + 1;
            iR = ir + 1;
            iC = ic + 1;
            plot(dX(iR,iC,iE),dY(iR,iC,iE),'o','color',colors(iE,:));
            data = [data;dX(iR,iC,iE),dY(iR,iC,iE)];
        end
    end
end
model = fitlm(data(:,1),data(:,2))
angle = atand(model.Coefficients.Estimate(2))
yOffSet = 500/6144*69000*model.Coefficients.Estimate(2)
