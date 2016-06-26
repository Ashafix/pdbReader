import numpy as np
import matplotlib.pyplot as plt
plt.xlim(0,20)
plt.ylim(0,1)
plt.ion()
y = []
i = 0
while i<10:
    temp = np.random.random()
    i += 1
    y.append(temp)
    if i>20:
        plt.xlim(i-20,i)
    plt.plot(y)
    plt.pause(0.1)
