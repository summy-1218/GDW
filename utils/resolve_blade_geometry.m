function [beamModel, r, c, theta, thickness, dr] = resolve_blade_geometry(params)
beamModel = struct();
if isfield(params,'beamModel')
    beamModel = params.beamModel;
elseif isfield(params,'beamModelFile')
    S = load(params.beamModelFile);
    if isfield(S,'beamModel')
        beamModel = S.beamModel;
    else
        beamModel = S;
    end
elseif exist('beamModel.mat','file')
    S = load('beamModel.mat');
    if isfield(S,'beamModel')
        beamModel = S.beamModel;
    else
        beamModel = S;
    end
end

% 读取几何参数：优先params，其次beamModel
if isfield(params,'r')
    r = params.r(:);
elseif isfield(beamModel,'distance')
    r = beamModel.distance(:);
else
    error('缺少叶片剖面径向位置: 请提供 params.r 或 beamModel.distance');
end

if isfield(params,'chord')
    c = params.chord(:);
elseif isfield(beamModel,'chord')
    c = beamModel.chord(:);
else
    error('缺少弦长分布: 请提供 params.chord 或 beamModel.chord');
end

if isfield(params,'twist')
    theta = params.twist(:);
elseif isfield(beamModel,'twist')
    theta = beamModel.twist(:);
else
    error('缺少扭角分布: 请提供 params.twist 或 beamModel.twist');
end

if isfield(params,'thickness')
    thickness = params.thickness(:);
elseif isfield(beamModel,'Thickness')
    thickness = beamModel.Thickness(:);
elseif isfield(beamModel,'thickness')
    thickness = beamModel.thickness(:);
else
    error('缺少翼型厚度参数: 请提供 params.thickness 或 beamModel.Thickness');
end

% 统一厚度单位(若为比例则转为百分数)
if max(thickness) <= 1.5
    thickness = thickness * 100;
end

% 叶素宽度由相邻距离决定
if isfield(params,'dr')
    dr = params.dr(:);
elseif isfield(beamModel,'distance')
    dr = [diff(r); diff(r(end-1:end))];
else
    dr = [diff(r); diff(r(end-1:end))];
end

% 匹配性检查
n = numel(r);
if numel(c) ~= n || numel(theta) ~= n || numel(thickness) ~= n || numel(dr) ~= n
    error('剖面数据长度不一致: r(%d), chord(%d), twist(%d), thickness(%d), dr(%d)', n, numel(c), numel(theta), numel(thickness), numel(dr));
end

% 若params与beamModel同时存在则校验一致性
if isfield(params,'r') && isfield(beamModel,'distance')
    if max(abs(params.r(:) - beamModel.distance(:))) > 1e-6
        error('params.r 与 beamModel.distance 不一致');
    end
end
if isfield(params,'chord') && isfield(beamModel,'chord')
    if max(abs(params.chord(:) - beamModel.chord(:))) > 1e-6
        error('params.chord 与 beamModel.chord 不一致');
    end
end
if isfield(params,'twist') && isfield(beamModel,'twist')
    if max(abs(params.twist(:) - beamModel.twist(:))) > 1e-6
        error('params.twist 与 beamModel.twist 不一致');
    end
end
if isfield(params,'thickness') && (isfield(beamModel,'Thickness') || isfield(beamModel,'thickness'))
    bmThick = thickness;
    if isfield(beamModel,'Thickness')
        bmThick = beamModel.Thickness(:);
    elseif isfield(beamModel,'thickness')
        bmThick = beamModel.thickness(:);
    end
    if max(bmThick) <= 1.5
        bmThick = bmThick * 100;
    end
    if max(abs(params.thickness(:) - bmThick(:))) > 1e-6
        error('params.thickness 与 beamModel 厚度数据不一致');
    end
end
end