import matplotlib.pyplot as plt
fig, ax = plt.subplots()

def arrow_(ax, plt, x, y, dx, dy, **kwargs):
    ax.arrow(x, y, dx, dy, **kwargs)
    plt.plot([x, x + dx + 0.1], [y, y + dx + 0.1], alpha=0)

arrow_(ax, plt, 1, 1, 1, 1)
ax.relim()
#plt.plot([1.8], [1.5], 'ro')
ax.autoscale_view()
plt.show()


