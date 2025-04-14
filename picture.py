import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# 网格参数
m, n = 20, 16

# 读取二进制h.grd数据
def read_fortran_unformatted(filename):
    data = []
    with open(filename, 'rb') as f:
        while True:
            header = np.fromfile(f, dtype=np.int32, count=1)
            if not header: break
            record = np.fromfile(f, dtype=np.float32, count=header[0]//4)
            footer = np.fromfile(f, dtype=np.int32, count=1)
            data.append(record)
    return data

# 读取位势高度场（取最后一个时间步）
hgt_data = read_fortran_unformatted('h.grd')
hgt = hgt_data[-1].reshape(n, m).T  # 转置为(m, n)

# 读取风场数据
def read_wind(filename):
    return np.loadtxt(filename).reshape(n, m).T  # 转置为(m, n)

u = read_wind('uc.dat')
v = read_wind('vc.dat')

# 清洗位势高度场数据
hgt = np.nan_to_num(hgt, nan=0.0)  # 替换NaN为0
hgt = np.clip(hgt, -1e10, 1e10)    # 限制极大/极小值

# 清洗风场数据
u = np.nan_to_num(u, nan=0.0)
v = np.nan_to_num(v, nan=0.0)
u = np.clip(u, -1000, 1000)
v = np.clip(v, -1000, 1000)
# 创建网格坐标
x = np.arange(m)
y = np.arange(n)
X, Y = np.meshgrid(x, y)

# 设置风场采样间隔（每隔2个点显示）
stride = 1  # 采样间隔
# 计算采样索引
skip = (slice(None, None, stride), slice(None, None, stride))

plt.figure(figsize=(12, 8))

# 创建自定义颜色映射: 蓝色-白色-橙色
colors = [(0.2, 0.4, 0.8), (1, 1, 1), (1, 0.6, 0.2)]  # 蓝色、白色、橙色（透明度降低）
cmap_name = 'blue_white_orange'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=256)

# 绘制位势高度场
cf = plt.contourf(X, Y, hgt.T,  # 转置回(n, m)方向
                 cmap=custom_cmap, 
                 levels=np.linspace(hgt.min(), hgt.max(), 20),
                 alpha=0.8)  # 降低一点不透明度，使颜色不太深
plt.colorbar(cf, label='Geopotential Height (m)')
# 添加等值线并每隔100标注数据
# 计算以100为间隔的等值线级别
min_level = np.floor(hgt.min() / 100) * 100
max_level = np.ceil(hgt.max() / 100) * 100
levels = np.arange(min_level, max_level + 100, 100)

# 绘制透明等值线
cs = plt.contour(X, Y, hgt.T,  # 转置回(n, m)方向
                 levels=levels, 
                 colors='black', 
                 linewidths=0.5,
                 alpha=1.0)
# 添加等值线标签
# 添加标注
plt.clabel(cs, cs.levels, inline=True, fmt='%d', fontsize=8)

# 绘制风场箭头
wind_plot = plt.quiver(X[skip], Y[skip], 
                      u.T[skip], v.T[skip],  # 转置后索引
                      scale=450,             # 调整箭头密度
                      color='black',
                      width=0.002,
                      headwidth=3,
                      headlength=4)

# 添加风场参考箭头
plt.quiverkey(wind_plot, 0.85, 0.92, 20, 
             label='20 m/s', 
             labelpos='E',
             color='black',
             coordinates='figure')

# 美化图形
plt.title('Geopotential Height with Wind Vectors')
plt.xlabel('X Grid Index')
plt.ylabel('Y Grid Index')
plt.grid(linestyle='--', alpha=0.5)

# 显示结果
plt.show()