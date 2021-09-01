optical文件夹和acoustic文件夹分别放置了光学仿真结果和声学仿真结果

【注】optical_skin_X.mat 中X表示血液模型中血糖的浓度，单位为mg/L（毫克每升）

optical_skin_X.mat表示光学仿真结果，内部变量说明如下：
（1） kgrid: 仿真用介质网格（正方形网格）--> 光学和声学的网格是共用的
（2） source: 初始压强分布
（3） absorb_energy: 吸收光子密度分布
（4） vmcmedium: 介质光学参数设置
         4.1 absorption_coefficient: 吸收系数；
         4.2 scattering_coefficient: 散射系数；
         4.3 scattering_anisotropy: 各向异性因子；
         4.4 refractive_index: 折射率；
（5） solution：光子密度分布

acoustic_skin_X.mat表示声学仿真结果，内部变量说明如下：
（1） source: 初始压强分布
（2） absorb_energy: 吸收光子密度分布
（3） sensor_data: 超声仿真结果（共21个探头，取最中间的列数据）
（4） vmcmedium: 介质光学参数设置
         4.1 absorption_coefficient: 吸收系数；
         4.2 scattering_coefficient: 散射系数；
         4.3 scattering_anisotropy: 各向异性因子；
         4.4 refractive_index: 折射率；
（5） solution：光子密度分布
（6） kgrid: 仿真用介质网格（正方形网格）--> 光学和声学的网格是共用的 （可能暂时还没有加进去）