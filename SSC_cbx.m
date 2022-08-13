%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Chen Bingxu, a 2020 graduate student of the School of Computer and Information Engineering, Henan University.
% 2022,08,13
%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;

tic
load('SLC_airport_20RFI.mat');
data_interference = S;
clear S;

Nr=8192;%距离向采样点数.totalBandwidth/Fr*Nr为有效带宽的采样点数
Fr=4.691840280000000e+07;%距离向采样率（过采样率
totalBandwidth=4.218944400000000e+07;%距离向处理总带宽

f_move = -0.93765e7;
t = ( -Nr/2 : (Nr/2-1) )/Fr;                % 距离时间轴
Sc = exp(1j*2*pi*f_move*t);
data_interference = data_interference.*Sc;

%距离向FFT,得到频谱
data_i_f = fty(data_interference);
clear data_interference;

%频谱切片图
fre_mean = mean(abs(data_i_f),1);
figure;
frequency_axis = linspace(-Fr/2,Fr/2,Nr);
plot(frequency_axis,fre_mean);
title('距离频率幅度谱');
xlabel('距离频域');
ylabel('幅度');
%% 选取截取的带宽频带
use_bandWidth = 4.043504928808594e+07;
rect = rectpuls(frequency_axis,use_bandWidth);%矩形窗过滤函数，过滤掉有效带宽之外的带宽,取4.043504928808594e+07
figure;
plot(fre_mean.*rect);

fre_mean_use = fre_mean(567:7626);%截取的有效采样点7060个点
use_axis = linspace(-use_bandWidth/2,use_bandWidth/2,7060);%截取的频率轴
figure;
plot(use_axis,fre_mean_use);
title('截取的距离频率幅度谱');
xlabel('距离频域');
ylabel('幅度');

data_i_f_use = data_i_f(:,567:7626);%截取的频谱
clear data_i_f;
%% 频谱校正
load('quxie_matrix.mat');%去斜矩阵

data_i_f_use_quxie = data_i_f_use.*quxie;%校正后的频谱
clear data_i_f_use quxie;
figure;
plot(use_axis,mean(abs(data_i_f_use_quxie),1));
title('频谱校正后的距离频率幅度谱');
xlabel('距离频域');
ylabel('幅度');
figure;
plot(mean(abs(data_i_f_use_quxie),1));
title('频谱校正后的距离频率幅度谱');
xlabel('距离频域');
ylabel('幅度');
%% 子带划分
% 20%
sub1_f = data_i_f_use_quxie(:,4702:6237);%干扰子带
sub2_f = data_i_f_use_quxie(:,3007:4542);%无干扰子带

sub1_f = [zeros(8192,4701),sub1_f,zeros(8192,823)];
sub2_f = [zeros(8192,3006),sub2_f,zeros(8192,2518)];

sub1_f = [zeros(8192,566),sub1_f,zeros(8192,566)];
sub2_f = [zeros(8192,566),sub2_f,zeros(8192,566)];
data_i_f_use_quxie = [zeros(8192,566),data_i_f_use_quxie,zeros(8192,566)];

%距离向IFFT
sub1 = ifty(sub1_f);
clear sub1_f;
sub2 = ifty(sub2_f);
clear sub2_f ;

data_i_use_quxie = ifty(data_i_f_use_quxie);
clear data_i_f_use_quxie;

%功率谱
S_sub1 = abs(sub1);
S_sub2 = abs(sub2);

S_data_interference = abs(data_i_use_quxie);
clear sub1 sub2 data_i_use_quxie;

%对消
S_rfi_1 = S_sub1 - S_sub2;

S_no_rfi = S_data_interference - S_rfi_1;

fenzi = S_data_interference-S_no_rfi;
fenzi_normF = norm(fenzi,'fro');
fenmu_normF = norm(S_data_interference,'fro');
RMSE_val = fenzi_normF/fenmu_normF

imwrite((uint16(65535*abs(S_no_rfi)/max(abs(S_no_rfi(:))))),'airport_noRFI.tiff');
toc