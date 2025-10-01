import numpy   as np            # Array operations.
import matplotlib.pyplot as plt # Plot data
import sys

for file in sys.argv[1:]:
    a = np.loadtxt(file)
    plt.plot(np.abs(a[:,-1]), a[:,0],label=file)

plt.legend()
plt.show()
