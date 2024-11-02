import pandas as pd
import matplotlib.pyplot as plt
import imageio
import numpy as np
import sys

def create_n_body_simulation_gif(csv_file_path, output_gif, grid_dim):
    # Load the CSV file
    try:
        data = pd.read_csv(csv_file_path)
    except Exception as e:
        print(f"Error reading the CSV file: {e}")
        return

    # Check and clean column names
    data.columns = data.columns.str.strip()

    # Verify that all required columns are present
    required_columns = ['iter_number', 'body_id', 'x_pos', 'y_pos', 'mass']
    if not all(col in data.columns for col in required_columns):
        print(f"Error: Missing one or more required columns: {', '.join(required_columns)}")
        return

    # Convert necessary columns to numeric types
    data['mass'] = pd.to_numeric(data['mass'], errors='coerce')
    data['x_pos'] = pd.to_numeric(data['x_pos'], errors='coerce')
    data['y_pos'] = pd.to_numeric(data['y_pos'], errors='coerce')

    # Extract unique iterations and body IDs
    iterations = data['iter_number'].unique()
    body_ids = data['body_id'].unique()

    # Preprocess data: create a dictionary with data for each iteration
    data_dict = {iter_num: iter_data for iter_num, iter_data in data.groupby('iter_number')}

    # Pre-calculate marker sizes
    marker_sizes = {body_id: mass for body_id, mass in zip(data['body_id'],data['mass'])}

    # Set up the plot
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(0, grid_dim)
    ax.set_ylim(0, grid_dim)
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
        ax.set_xlim(0, grid_dim)
        ax.set_ylim(0, grid_dim)
        ax.set_aspect('equal')
        ax.set_title(f'N-Body Simulation - time elapsed {iter_num/4:.2f}s')
        ax.set_xlabel('X Position')
        ax.set_ylabel('Y Position')

        # Get the data for the current iteration
        iter_data = data_dict[iter_num]
        if iter_data.empty:
            continue  # Skip this iteration if there's no data

        # Plot each body
        for body_id in body_ids:
            if body_id in iter_data['body_id'].values:  # Check if body_id exists in the current iteration
                body_data = iter_data[iter_data['body_id'] == body_id].iloc[0]
                x_pos = body_data['x_pos']
                y_pos = body_data['y_pos']

                # Set marker size, defaulting to 5 if the body_id is not in marker_sizes
                marker_size = marker_sizes.get(body_id, 5);
                
                # Plot the body
                ax.plot(x_pos, y_pos, 'o', markersize=marker_size, label=f'Body {body_id}' if iter_num == iterations[0] else "")

        # Save frame
        fig.canvas.draw()
        image = np.array(fig.canvas.buffer_rgba())
        frames.append(image)

    # Create GIF
    imageio.mimsave(output_gif, frames, fps=10, quality=8)
    plt.close(fig)

    print(f"GIF saved as {output_gif}")

# If you want to run the function from the command line:
if __name__ == "__main__":
    csv_file_path = sys.argv[1] if len(sys.argv) > 1 else 'output_simulation.csv'
    gif_file_path = sys.argv[2] if len(sys.argv) > 2 else 'output_simulation.gif'
    grid_dim = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    create_n_body_simulation_gif(csv_file_path, gif_file_path, grid_dim)
