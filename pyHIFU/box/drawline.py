import numpy as np
import matplotlib.pyplot as plt
# import seaborn
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d
from random import random
import matplotlib.ticker as plticker


def plot_straight_line(se, points, ax):
    p1 = se[0]
    p2 = se[1]
    # t0 = time.time()
    # print("optimized:", time.time() - t0)

    d = points.max() - points.min() + 1
    xmin = points[:,0].min()
    ymin = points[:,1].min()

    data = np.zeros((d, d, 3))
    for p in points:
        data[p[0]-xmin, p[1]-ymin, :] = [1,1,1]
    ax.imshow(data)
    ax.plot([p1[1]-ymin, p2[1]-ymin], [p1[0]-xmin, p2[0]-xmin])
    ax.scatter(points[:,1]-ymin, points[:,0]-xmin, s=1, color="red")

    ticks = np.arange(0, d, 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticks+xmin)
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks+ymin)
    ax.grid(b=True, which='both')
    

def bresenham(p1, p2):
    """
    Bresenham line drawing algorithm
    `p1`, `p2`: np.array()
    return: list of all the points to plot
    Assuming: slope `m` in [0, 1]
    """
    dim = 2

    x1f = p1[0]
    x2f = p2[0]

    y1f = p1[1]
    y2f = p2[1]
    
    slope = (y2f - y1f) / (x2f - x1f)
    x1 = round(x1f)
    y1 = round(y1f)
    x2 = round(x2f)

    points = list()

    # initial error
    err = (slope*(x1 - x1f) + y1f) - y1
    y = y1
    for x in range(x1, x2+1):
        points.append([x, y])
        err_ = err+slope
        # print(x,err_ < 0.5, err, err_)
        if err_ < 0.5:
            y = y
            err = err_
        else:
            # also select the pixel neglected at turning points
            err_mid = err + 0.5*slope
            if err_mid < 0.5:
                points.append([x+1, y])
            else:
                points.append([x, y+1])

            y = y + 1
            err = err_ - 1
    return np.array(points)

def bresenham3D(p1, p2):
    pass

def bresenham2D(p1, p2):
    # Optimized to use integer addition only
    # though in python it seems float is faster than int
    t0 = time.time()
    dim = 2

    x1f = p1[0]
    x2f = p2[0]
    y1f = p1[1]
    y2f = p2[1]
    
    # all np.integers from below
    x1 = np.int(x1f)
    y1 = np.int(y1f)
    x2 = np.int(x2f)
    y2 = np.int(y2f)
    
    dx = x2 - x1
    dy = y2 - y1
    m = 2*dy
    mx = 2*dx

    points = np.ndarray((x2-x1+1, dim), dtype=np.int)
    e = 0
    y = y1
    for x in range(x1, x2+1):
        points[x-x1,:] = [x, y]
        e = e + m
        if e < dx:
            y = y
        else:
            y = y + 1
            e = e - mx
    print("optimized:", time.time() - t0)
    # return points
    

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111)
    p1 = [-7.6, -2.1]
    p2 = [8987654.3, 4986.9]
    points = bresenham(p1, p2)
    plot_straight_line([p1, p2], points, ax)
    plt.show()