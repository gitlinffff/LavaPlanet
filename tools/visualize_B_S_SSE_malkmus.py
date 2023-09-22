import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

# Load data from CSV
data = pd.read_csv('E:/IceGiantProject/B1_BS_SSE.csv')  # Replace with your actual file name

# Extract x, y, z columns from the data
x = data['B']
y = data['S']
z = data['SSE']

# Create a figure
fig = plt.figure()
ax = fig.add_subplot()

# Plot the data
ax.scatter(x, y, c=z, cmap='hot',marker='.')  # Adjust colors and markers as needed

ax.set_xscale('log')
ax.set_yscale('log')

# Set labels
ax.set_xlabel('B')
ax.set_ylabel('S')
#ax.set_zlabel('Z Label')

# Show the plot
plt.show()
#output_fig = '/home/linfel/LavaPlanet/images/B1_3d.png'
#plt.savefig(output_fig,dpi=300)
