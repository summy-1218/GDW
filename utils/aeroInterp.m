function [cl,cd,cm] = aeroInterp(thicknessi,alphai)
persistent DU21 DU25 DU30 DU35 DU40 DU90 cachedAlpha cachedX cachedY ...
           z_cl z_cd z_cm;
if isempty(DU21)
    load cl_cd_cm.mat DU21 DU25 DU30 DU35 DU40 DU90;
    alpha_ref = DU21(:,1);
    % 厚度轴以比例(0~1)存储; 调用方须以百分比(%)传入 thicknessi, 函数内除以100转换
    % 例: thicknessi=21 对应 t/c=0.21
    thickness_pts = [0.1,0.21,0.25,0.3,0.35,0.4,0.9,1.0];
    cl_21 = pchip(DU21(:,1),DU21(:,2),alpha_ref); 
    cd_21 = pchip(DU21(:,1),DU21(:,3),alpha_ref);
    cm_21 = pchip(DU21(:,1),DU21(:,4),alpha_ref);

    cl_10 = cl_21;
    cd_10 = cd_21;
    cm_10 = cm_21;

    cl_25 = pchip(DU25(:,1),DU25(:,2),alpha_ref);
    cd_25 = pchip(DU25(:,1),DU25(:,3),alpha_ref);
    cm_25 = pchip(DU25(:,1),DU25(:,4),alpha_ref);
    
    cl_30 = pchip(DU30(:,1),DU30(:,2),alpha_ref);
    cd_30 = pchip(DU30(:,1),DU30(:,3),alpha_ref);
    cm_30 = pchip(DU30(:,1),DU30(:,4),alpha_ref);
    
    cl_35 = pchip(DU35(:,1),DU35(:,2),alpha_ref);
    cd_35 = pchip(DU35(:,1),DU35(:,3),alpha_ref);
    cm_35 = pchip(DU35(:,1),DU35(:,4),alpha_ref);
    
    cl_40 = pchip(DU40(:,1),DU40(:,2),alpha_ref);
    cd_40 = pchip(DU40(:,1),DU40(:,3),alpha_ref);
    cm_40 = pchip(DU40(:,1),DU40(:,4),alpha_ref);
    
    cl_90 = pchip(DU90(:,1),DU90(:,2),alpha_ref);
    cd_90 = pchip(DU90(:,1),DU90(:,3),alpha_ref);
    cm_90 = pchip(DU90(:,1),DU90(:,4),alpha_ref);
    
    cl_100 = cl_90;
    cd_100 = cd_90;
    cm_100 = cm_90;

    [cachedX, cachedY] = meshgrid(thickness_pts, alpha_ref);  % 转换为网格矩阵
    z_cl = [cl_10,cl_21,cl_25,cl_30,cl_35,cl_40,cl_90,cl_100];
    z_cd = [cd_10,cd_21,cd_25,cd_30,cd_35,cd_40,cd_90,cd_100];
    z_cm = [cm_10,cm_21,cm_25,cm_30,cm_35,cm_40,cm_90,cm_100];
    
    cachedAlpha = alpha_ref;
end
cl = interp2(cachedX,cachedY,z_cl,thicknessi/100,rad2deg(alphai),'linear');
cd = interp2(cachedX,cachedY,z_cd,thicknessi/100,rad2deg(alphai),'linear');
cm = interp2(cachedX,cachedY,z_cm,thicknessi/100,rad2deg(alphai),'linear');

end
