from trident import Trident
import matplotlib.pyplot as plt
import numpy as np
import cProfile

tr = Trident([1,1,1],[2,3,4],
            [1,1,1],[5,2,3],
            [1,1,1],[5,3,4])

t = np.arange(0.5, 10, 0.1)
area_list = np.ndarray(t.shape)

cProfile.run('''
for i,d in enumerate(t):
    print(i)
    area_list[i]=tr.get_area_at(d)
''')


plt.plot(np.square(t), area_list)
plt.show()