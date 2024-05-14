% clear;
% clc;

tic;
%折射率：587.6nm单波长下
sk2=1.6073788586;
sk16=1.6204079331;
f5=1.6034171782;
%折射率：555.0nm单波长下
bk7 = 1.5185223876;

wlUM = 587.6e-6;
air=1;
opt = opt_model();
a = opt.seq_model;
a.source_thick = 300; %inf=1e10
a.source_field_objheight = [5,30];
a.setEntrancePupil(33.33); %设置入瞳直径

%% 双高斯透镜（单波长587.6nm）
a.initialInterface()
a.addInterface('Conic', 54.153, 8.747, sk2, 1);
a.addInterface('Conic', 152.522, 0.5, air, 1);
a.addInterface('Conic', 35.951, 14, sk16, 1);
a.addInterface('Conic', inf, 3.777, f5, 1);
a.addInterface('Conic', 22.270, 14.253, air, 1);
a.addInterface('Conic', inf, 12.428, air, 1);
a.set_stop();
a.addInterface('Conic', -25.685, 3.777, f5, 1);
a.addInterface('Conic', inf, 10.834, sk16, 1);
a.addInterface('Conic', -36.980, 0.5, air, 1);
a.addInterface('Conic', 196.417, 6.858, sk16, 1);
a.addInterface('Conic', -67.148, 94.071, air, 1);

a.field_height_x = a.source_field_objheight(1);
a.field_height_y = a.source_field_objheight(2);
a.update();

exp_pos = a.opt_model.ParaxialModel.ExitPupil_position;
exp_dia = a.opt_model.ParaxialModel.ExitPupil_diameter;
ray_in_enp = a.ray_in_enp./(a.opt_model.ParaxialModel.EntrancePupil_diameter/2);
s_in_img_plane = a.s_back;
s_in_img_plane(:,3) = -exp_pos;
optical_path_obj_and_img = a.optical_path_total;

single_ray = a.specific_trace(0,1); %子午面(YOZ)
cos1 = single_ray(6)/sqrt(single_ray(5)^2+single_ray(6)^2);
sin1 = sqrt(1 - cos1^2);
tfno = 1/(2*air*sin1); %子午工作F数

single_ray = a.specific_trace(1,0); %弧矢面(XOZ)
cos1 = single_ray(6)/sqrt(single_ray(5)^2+single_ray(6)^2);
sin1 = sqrt(1 - cos1^2);
sfno = 1/(2*air*sin1); %弧矢工作F数

wfno = min(sfno,tfno); %工作F数

% 反向追迹到出瞳
a.addInterface('Conic', inf, exp_pos, air, 1); %像面
a.addInterface('Conic', inf, 0, air, 1); %出瞳面
a.addInterface('Conic', -exp_pos, -exp_pos, air, 1); %参考球面

a.field_height_x = a.source_field_objheight(1);
a.field_height_y = a.source_field_objheight(2);
a.update();

%% 计算波前图
s_in_ideal_plane = a.opt_model.interfaces(end-1).s_inter;
optical_path_exp_and_img = (exp_pos/abs(exp_pos))*sqrt((s_in_ideal_plane(:,1) - s_in_img_plane(:,1)).^2 + (s_in_ideal_plane(:,2) - s_in_img_plane(:,2)).^2 + (s_in_ideal_plane(:,3) - s_in_img_plane(:,3)).^2);
optical_path_total = optical_path_obj_and_img + optical_path_exp_and_img;
OPD = optical_path_total - optical_path_total(1);
WF = -OPD./(0.5876e-3);

sampling = 128;
pupil_sample = sampling - 1;
x_detect_spot = linspace(-1,1,pupil_sample);
y_detect_spot = linspace(-1,1,pupil_sample);
wfmap_color = zeros(pupil_sample,pupil_sample);
pixel_ray_count = zeros(pupil_sample,pupil_sample);
for j = 1:size(ray_in_enp,1)
    [~,x_near] = min(abs(x_detect_spot(:) - ray_in_enp(j,1)));
    [~,y_near] = min(abs(y_detect_spot(:) - ray_in_enp(j,2)));
    wfmap_color(y_near,x_near)  = wfmap_color(y_near,x_near) + WF(j);
    pixel_ray_count(y_near,x_near)  = pixel_ray_count(y_near,x_near) + 1;
end
wfmap_color = wfmap_color./pixel_ray_count;

PV = max(wfmap_color,[],'all')-min(wfmap_color,[],'all');
wave_mean = mean(wfmap_color,[1 2],"omitnan");
RMS = rms((wfmap_color-wave_mean),[1 2],"omitnan");

wavefront_map = wfmap_color;
% wavefront_map = flipud(wavefront_map);
wf_x_label = linspace(-1, 1, 3);
wf_y_label = linspace(-1, 1, 3);

figure;
imagesc(wf_x_label,wf_y_label,wavefront_map);
axis equal tight;
axis xy;

h1 = colorbar;
t = get(h1,'Limits');
T1 = linspace(t(1),t(2),11);
set(h1,'Ticks',T1);
TL1 = arrayfun(@(x) sprintf('%.2f',x),T1,'un',0);
set(h1,'TickLabels',TL1);

load falsecolor.mat;
colormap(CustomColormap);
xticks(wf_x_label);
yticks(wf_y_label);
title({'Calculated Wavefront Map';['Wavelehgth = ',num2str(wlUM*1e6),'nm, field = [',num2str(a.source_field_objheight(1)),', ',num2str(a.source_field_objheight(2)),']mm']});
xlabel('X-光瞳（相对单位）');
ylabel('Y-光瞳（相对单位）');
set(get(h1,'Title'),'string','Waves');
disp('Wavefront map: ');
disp(['  PV = ', num2str(PV),' waves']);
disp(['  RMS = ', num2str(RMS),' waves']);
disp(['  Exit pupil diameter = ', num2str(exp_dia),' mm']);

%% 波前到PSF
% 将光程差转换为波
waves = 2*pi*wfmap_color;
real = cos(waves);
imag = sin(waves);
real(isnan(real)) = 0;
imag(isnan(imag)) = 0;
% 填充0
[ns,~] = size(wfmap_color);
padding = ns*2;
real = padarray(real,[padding padding],0,'both');
imag = padarray(imag,[padding padding],0,'both');
% 计算复数光瞳函数
pf = real + 1i*imag;
% 得到FFT复振幅
amp = ifftshift(fft2(fftshift(pf)));
% 得到强度
psf = abs(amp.^2);
% 缩放到Strehl Ratio
numPassed = nnz(~isnan(wfmap_color));
sr = 1 / (numPassed * numPassed);
psf = psf * sr;
% 旋转方向
psf = flipud(psf);
psf = fliplr(psf);
% PSF绘制范围
[px, py] = size(psf);
factor = 1;
zoom = 1 / sqrt(32 / ns);
psf_zoom = psf(fix(px/2) - fix(px/(zoom*2)):fix(px/2) + fix(px/(zoom*2)), fix(py/2) - fix(py/(zoom*2)):fix(py/2) + fix(py/(zoom*2)))*factor;
delta = wlUM * wfno * ((sampling-2) * (sqrt(32/sampling) / (sampling*2)));
% PSF绘制
img_size = 2*ns*delta*1e3; %单位：微米
psf_x_label = linspace(-ns*delta, ns*delta, 3)*1e3;
psf_y_label = linspace(-ns*delta, ns*delta, 3)*1e3;
figure;
imagesc(psf_x_label,psf_y_label,psf_zoom);
axis equal tight;
axis xy;

h2 = colorbar;
t = get(h2,'Limits');
T2 = linspace(t(1),t(2),11);
set(h2,'Ticks',T2);
TL2 = arrayfun(@(x) sprintf('%.4f',x),T2,'un',0);
set(h2,'TickLabels',TL2);

colormap(CustomColormap);
xticks(psf_x_label);
yticks(psf_y_label);
title({'Calculated FFT PSF';['Wavelehgth = ',num2str(wlUM*1e6),'nm, field = [',num2str(a.source_field_objheight(1)),', ',num2str(a.source_field_objheight(2)),']mm']});
xlabel('X位置: μm');
ylabel('Y位置: μm');
disp('FFT PSF: ');
disp(['  Side = ', num2str(img_size),' μm']);
disp(['  Reference Coordinates = ', num2str(s_in_img_plane(1,1)),', ',num2str(s_in_img_plane(1,2)),' mm']);
disp(['Strehl Ratio  : ' ,num2str(max(psf_zoom,[],'all'))]);
