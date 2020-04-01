## 如何使用SRW：核心是函数库的安装，及编程方式

### 1. 直接按照SRW python

   <https://srwpy.github.io/srwpy/index.html>
   <https://github.com/srwpy/srwpy>

    安装：pip install srwpy

    *使用*：
     >from srwpy.srwlib import *
     >from srwpy.uti_plot import *


### 2. 下载依赖库文件

   <https://github.com/ochubar/SRW>

    参照本文件夹下的SRW_SE_example

    *使用*：
    > from srwlib import *
    > from uti_plot import *


### 3. 安装OASYS，基于OASYS_SRW生成代码，运行程序

<https://www.aps.anl.gov/Science/Scientific-Software/OASYS>
   
    *使用*：
    > from oasys_srw.srwlib import *
    > from oasys_srw.uti_plot import *

### 4. 使用Sirepo

 <https://www.bnl.gov/ps/groups/rd/simulations/tools.php>

 <https://www.sirepo.com/light#/home>

 <https://github.com/radiasoft/sirepo>

    *使用*：
    > import srwl_bl
    > import srwlib
    > import srwlpy
