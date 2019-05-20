import matplotlib.pyplot as plt
import numpy as np
from pyHIFU.io.config import readjson
from matplotlib.ticker import FormatStrFormatter


n_rays_data = readjson(json_path="data/d_v_n_rays.json")
theta_max_data = readjson(json_path="data/d_v_theta_max.json")
trident_angle_data = readjson(json_path="data/d_v_trident_angle.json")

n_rays_l = list(n_rays_data.keys())
n_rays_l = np.array(n_rays_l, dtype=np.float64)
n_rays_v = n_rays_data.values()

theta_max_l = list(theta_max_data.keys())
theta_max_l = np.array(theta_max_l, dtype=np.float64)
theta_max_v = theta_max_data.values()

trident_angle_l = list(trident_angle_data.keys())
trident_angle_l = np.array(trident_angle_l, dtype=np.float64)
trident_angle_v = trident_angle_data.values()

fig, ax = plt.subplots()

# ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax.plot(n_rays_l, n_rays_v)
# ax.title.set_text("number of rays")
# ax.set_xlabel("n_rays")
# ax.set_ylabel("pnorm distance")

# ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax.plot(theta_max_l, theta_max_v)
# ax.title.set_text(r"$\theta_{max}$")
# ax.set_xlabel(r"$\theta_{max}$")
# ax.set_ylabel("pnorm distance")

ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.plot(trident_angle_l, trident_angle_v)
ax.title.set_text("trident angle")
ax.set_xlabel("trident angle")
ax.set_ylabel("pnorm distance")

plt.show()
