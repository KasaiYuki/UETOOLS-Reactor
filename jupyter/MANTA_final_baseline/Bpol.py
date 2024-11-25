import numpy as np
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Read data from text files
# Assuming your x, y coordinates are in 'coordinates.txt' and magnitudes have 3 columns in 'magnitudes.txt'

# Load x and y coordinates (assuming each value is space-separated)
coords = np.loadtxt('v9_b_grid.txt')

# Load magnitudes (assuming there are 3 columns, and you want to use only the first column)
magnitudes = np.loadtxt('v9_b_mag.txt')

# Step 2: Extract x, y, and magnitude columns
x = coords[:, 0]  # First column is x
y = coords[:, 1]  # Second column is y

# Extract only the first column of magnitudes
Br = np.array(magnitudes[:, 0])
Bz = np.array(magnitudes[:, 2])
Bpol = np.sqrt(Br**2+Bz**2)
print(Bpol)

# Step 3: Create a grid for plotting
# Create a grid of x and y values
grid_x, grid_y = np.meshgrid(np.unique(x), np.unique(y))

# Interpolate the magnitude values onto the grid
grid_magnitude = np.zeros_like(grid_x)

# Using nearest neighbor interpolation (you could use other methods like linear if needed)
for i in range(len(x)):
    xi = np.where(np.unique(x) == x[i])[0][0]
    yi = np.where(np.unique(y) == y[i])[0][0]
    grid_magnitude[yi, xi] = Bpol[i]
print(max(Bpol))
# Step 4: Create the color plot
plt.figure(figsize=(8, 6))
plt.scatter(x,y, c=Bpol)  # You can choose other color maps
plt.colorbar(label='Magnitude (First Column)')
plt.xlabel('X Coordinate')
plt.ylim(-2,2)
plt.xlim(1,3)
plt.clim(0,2)
plt.ylabel('Y Coordinate')
plt.title('2D Color Plot of Magnitudes (First Column)')

# Show the plot
plt.show()
