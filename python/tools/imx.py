import numpy as np
import matplotlib.pyplot as plt

class IndexTracker(object):
    def __init__(self, ax, X, cmap='gray'):
        self.ax = ax
        ax.set_title('use arrow keys to navigate')

        self.X = X
        rows, cols, self.slices = X.shape
        self.ind = self.slices//2

        self.cmap = cmap

        self.im = ax.imshow(self.X[:, :, self.ind], cmap=cmap)
        self.update()

    def onscroll(self, event):
        if event.key == 'up' or event.key == 'right':
            self.ind = (self.ind + 1) % self.slices
        elif event.key == 'left' or event.key == 'down':
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(self.X[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()

def imx(images, fig=None, **kwargs):
    if isinstance(images, list):
        N = len(images)

        if fig is None:
            fig, ax = plt.subplots(1, N)
        else:
            ax = fig.axes

        tracker = []
        for i, im in enumerate(images):
            tracker.append(IndexTracker(ax[i], im, cmap=kwargs.get('cmap', 'gray')))
            fig.canvas.mpl_connect('key_press_event', tracker[i].onscroll)
    else:
        if fig is None:
            fig, ax = plt.subplots(1, 1)
        else:
            ax = fig.axes[0]
        tracker = IndexTracker(ax, images, cmap=kwargs.get('cmap', 'gray'))

        fig.canvas.mpl_connect('key_press_event', tracker.onscroll)
    plt.show()

    # We must return these to avoid garbage collection
    return fig, tracker
