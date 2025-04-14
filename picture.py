import numpy as np
import matplotlib.pyplot as plt

# 网格参数
m, n = 20, 16

# 读取二进制数据文件
def read_fortran_unformatted(filename, dtype=np.float32):
    data = []
    with open(filename, 'rb') as f:
        while True:
            # 读取记录头（4字节长度）
            header = np.fromfile(f, dtype=np.int32, count=1)
            if not header:
                break
            # 计算元素个数
            n_elements = header[0] // np.dtype(dtype).itemsize
            # 读取数据
            record = np.fromfile(f, dtype=dtype, count=n_elements)
            # 读取记录尾（4字节长度）
            footer = np.fromfile(f, dtype=np.int32, count=1)
            data.append(record)
    return data

# 读取h.grd文件
data = read_fortran_unformatted('h.grd')

# 解析位势高度场数据（假设文件中只有一个时间层）
if len(data) == 1:
    hgt = data[0].reshape((n, m)).T  # Fortran列优先存储，转置为(m, n)
elif len(data) > 1:
    hgt = data[-1].reshape((n, m)).T  # 取最后一个时间步

# 创建网格坐标
x = np.arange(m)
y = np.arange(n)
X, Y = np.meshgrid(x, y)

# 绘制位势高度场
plt.figure(figsize=(10, 6))
plt.contourf(X, Y, hgt.T, cmap='jet', levels=20)  # 转置回(n, m)方向
plt.colorbar(label='Geopotential Height (m)')
plt.title('Forecast Geopotential Height Field')
plt.xlabel('X Grid Index')
plt.ylabel('Y Grid Index')
plt.grid(linestyle='--', alpha=0.5)
plt.show()