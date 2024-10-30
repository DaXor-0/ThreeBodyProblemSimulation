import pandas as pd
import matplotlib.pyplot as plt
import imageio
import numpy as np

# Load the CSV file
file_path = 'output_simulation.csv'
data = pd.read_csv(file_path)

# Clean column names
data.columns = data.columns.str.strip()

# Define grid limits
GRID_MIN, GRID_MAX = 0, 100

# Extract unique iterations
iterations = data['iter_number'].unique()

# Group data by iteration and body ID for efficient access
grouped = data.groupby(['iter_number', 'body_id']).agg({'x_pos': 'first', 'y_pos': 'first', 'mass': 'first'}).reset_index()

# Set up the plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(GRID_MIN, GRID_MAX)
ax.set_ylim(GRID_MIN, GRID_MAX)
ax.set_aspect('equal')
ax.set_title('N-Body Simulation')
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')

# Prepare to save frames as images
frames = []

# Iterate over each time step
for iter_num in iterations:
    # Filter grouped data for the current iteration
    iter_data = grouped[grouped['iter_number'] == iter_num]

    # Plot the bodies
    ax.clear()  # Clear only the contents of the plot
    ax.set_xlim(GRID_MIN, GRID_MAX)
    ax.set_ylim(GRID_MIN, GRID_MAX)
    ax.set_title(f'N-Body Simulation - Iteration {iter_num}')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')

    # Scatter plot for all bodies
    ax.scatter(iter_data['x_pos'], iter_data['y_pos'], s=5, alpha=0.6)

    # Save frame
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8').reshape(fig.canvas.get_width_height()[::-1] + (3,))
    frames.append(image)

# Create GIF
output_gif = 'n_body_simulation.gif'
imageio.mimsave(output_gif, frames, fps=10)
plt.close(fig)

print(f"GIF saved as {output_gif}")
