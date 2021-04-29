# sound-source-localization-algorithm_DOA_estimation
* 语音信号处理的宽带说话人(声源)定位（DOA估计）算法

**Abstract** 本仓库是面向语音信号的声源定位传统算法

**关键词**：声源定位(sound source localization)、DOA估计（DOA estimation）、TDOA估计（TDOA estimation）、麦克风阵列信号处理(microphone array signal processing)
## ssl_tools
包含SRP-PHAT(GCC-PHAT)、MUSIC、beamforming（波束形成）三类算法
*  SRP：SRP-PHAT、非线性SRP-PHAT
*  MUSIC
*  beamforming:基于延迟求和(DS)的SNR方位谱估计、基于MVDR的SNR方位谱估计及其对应的频率加权改进算法


## 与语音信号处理的宽带声源定位相关的参考资源
### 竞赛
* acoustic source LOCalization And TrAcking [[LOCATA]](https://locata.lms.tf.fau.de/)
* Detection and Classification of Acoustic Scenes and Events [[DCASE]](http://dcase.community/challenge2020/task-sound-event-localization-and-detection)
### 多通道数据集生成算法
* rir-generator [[Code]](https://github.com/ehabets/RIR-Generator)
* ROOMSIM[[Code]](https://github.com/Wenzhe-Liu/ROOMSIM)
### 开源代码 
#### 基于时延的定位
* A simple DOA GUI 
[[Code]](https://github.com/wangwei2009/DOA)
#### 基于波束形成的定位
* DNN_Localization_And_Separation 
[[Code]](https://github.com/shaharhoch/DNN_Localization_And_Separation)
#### 双耳定位
* binauralLocalization 
[[Code]](https://github.com/nicolasobin/binauralLocalization)
* Binaural-Auditory-Localization-System 
[[Code]](https://github.com/r04942117/Binaural-Auditory-Localization-System)
* Binaural_Localization:ITD-based localization of sound sources in complex acoustic environments [[Code]](https://github.com/Hardcorehobel/Binaural_Localization)
#### 高分辨率定位
* WSCM-MUSIC
[[Code]](https://github.com/xuchenglin28/WSCM-MUSIC)
#### 基于聚类定位
* messl:Model-based EM Source Separation and Localization 
[[Code]](https://github.com/mim/messl) [[Paper]](https://www.ee.columbia.edu/~ronw/pubs/taslp09-messl.pdf) 
* fast_sound_source_localization_using_TLSSC:Fast Sound Source Localization Using Two-Level Search Space Clustering
[[Code]](https://github.com/LeeTaewoo/fast_sound_source_localization_using_TLSSC)
#### 窄带定位
* doa-tools
[[Code]](https://github.com/morriswmz/doa-tools)
* 麦克风声源定位 [[Code]](https://github.com/xiaoli1368/Microphone-sound-source-localization)

