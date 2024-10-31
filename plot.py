import pandas as pd
import matplotlib.pyplot as plt
import imageio
import numpy as np
import sys

def create_n_body_simulation_gif(csv_file_path, output_gif='n_body_simulation.gif', grid_max=100):
    # Load the CSV file
    data = pd.read_csv(csv_file_path)

    # Check and clean column names
    data.columns = data.columns.str.strip()  # Remove extra whitespace around column names if present

    # Define grid limits
    GRID_MIN = 0

    # Extract unique iterations and body IDs
    iterations = data['iter_number'].unique()
    body_ids = data['body_id'].unique()

    # Set up the plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(GRID_MIN, grid_max)
    ax.set_ylim(GRID_MIN, grid_max)
    ax.set_aspect('equal')
    ax.set_title('N-Body Simulation')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')

    # Prepare to save frames as images
    frames = []

    # Iterate over each time step
    for iter_num in iterations:
        # Clear plot for the current frame
        ax.cla()
        ax.set_xlim(GRID_MIN, grid_max)
        ax.set_ylim(GRID_MIN, grid_max)
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
                marker_size = mass * 10

                # Plot the body
                ax.plot(x_pos, y_pos, 'o', markersize=marker_size, label=f'Body {body_id}' if iter_num == 0 else "")

        # Save frame
        fig.canvas.draw()
        image = np.array(fig.canvas.buffer_rgba())
        frames.append(image)

    # Create GIF
    imageio.mimsave(output_gif, frames, fps=15)
    plt.close(fig)

    print(f"GIF saved as {output_gif}")

# If you want to run the function from the command line:
if __name__ == "__main__":
    csv_file_path = sys.argv[1] if len(sys.argv) > 1 else 'output_simulation.csv'
    grid_max = int(sys.argv[3]) if len(sys.argv) > 2 else 100
    create_n_body_simulation_gif(csv_file_path, grid_max=grid_max)
