import pandas as pd
import matplotlib.pyplot as plt
import imageio
import numpy as np

# Load the CSV file
file_path = 'output_simulation.csv'
data = pd.read_csv(file_path)

# Check and clean column names
print(data.columns)  # Debugging step
data.columns = data.columns.str.strip()  # Remove extra whitespace around column names if present

# Define grid limits based on provided GRID_MIN and GRID_MAX
GRID_MIN, GRID_MAX = 0, 100

# Extract unique iterations and body IDs
iterations = data['iter_number'].unique()
body_ids = data['body_id'].unique()

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
# Iterate over each time step
for iter_num in iterations:
    # Clear plot for the current frame
    ax.cla()
    ax.set_xlim(GRID_MIN, GRID_MAX)
    ax.set_ylim(GRID_MIN, GRID_MAX)
    ax.set_aspect('equal')
    ax.set_title(f'N-Body Simulation - Iteration {iter_num}')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')

    # Filter data for the current iteration
    iter_data = data[data['iter_number'] == iter_num]

    # Plot each body
    for body_id in body_ids:
        body_data = iter_data.loc[iter_data['body_id'] == body_id]  # Use .loc to filter
        if not body_data.empty:
            # Get position of the body
            x_pos = body_data['x_pos'].values[0]
            y_pos = body_data['y_pos'].values[0]
            mass = body_data['mass'].values[0]

            # Set size of the marker relative to mass
            marker_size = 5

            # Plot the body
            ax.plot(x_pos, y_pos, 'o', markersize=marker_size, label=f'Body {body_id}' if iter_num == 0 else "")

    # Save frame
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8').reshape(fig.canvas.get_width_height()[::-1] + (3,))
    frames.append(image)

# Create GIF
output_gif = 'n_body_simulation.gif'
imageio.mimsave(output_gif, frames, fps=10)
plt.close(fig)

print(f"GIF saved as {output_gif}")
