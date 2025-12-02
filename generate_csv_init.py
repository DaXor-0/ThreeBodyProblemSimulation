import csv
import random
import argparse

def generate_initial_conditions(
    num_bodies, 
    pos_range, 
    vel_range, 
    acc_range, 
    mass_range,
    zero_acc, 
    output_file
):
    """
    Generate initial conditions for an N-body simulation and save to a CSV file.

    Args:
        num_bodies (int): Number of bodies to simulate.
        pos_range (tuple): Range for initial positions (min, max).
        vel_range (tuple): Range for initial velocities (min, max).
        acc_range (tuple): Range for initial accelerations (min, max).
        mass_range (tuple): Range for masses (min, max).
        zero_acc (bool): If True, set initial accelerations to zero.
        output_file (str): File name for the output CSV.
    """
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write the header row to the CSV file
        writer.writerow(["x", "y", "v_x", "v_y", "a_x", "a_y", "mass"])

        for _ in range(num_bodies):
            # Generate random values for the x and y positions within the specified range
            x = random.uniform(*pos_range)
            y = random.uniform(*pos_range)

            # Generate random values for the x and y components of velocity within the specified range
            v_x = random.uniform(*vel_range)
            v_y = random.uniform(*vel_range)

            # If zero_acc is True, set accelerations to zero; otherwise, generate random values for acceleration
            if zero_acc:
                a_x, a_y = 0.0, 0.0
            else:
                a_x = random.uniform(*acc_range)
                a_y = random.uniform(*acc_range)

            # Generate a random value for the mass within the specified range
            mass = random.uniform(*mass_range)

            # Write the generated data to the CSV file as a new row
            writer.writerow([x, y, v_x, v_y, a_x, a_y, mass])

    # Inform the user that the file has been successfully created
    print(f"Initial conditions saved to {output_file}")

if __name__ == "__main__":
    # Create an argument parser to handle command-line arguments
    parser = argparse.ArgumentParser(description="Generate initial conditions for an N-body simulation.")

    # Add arguments to specify the number of bodies and the ranges for positions, velocities, and accelerations
    parser.add_argument("-n", "--num_bodies", type=int, required=True, help="Number of bodies")
    parser.add_argument("-p", "--pos_range", type=float, nargs=2, required=True, help="Range for positions (min max)")
    parser.add_argument("-v", "--vel_range", type=float, nargs=2, required=True, help="Range for velocities (min max)")
    parser.add_argument("-a", "--acc_range", type=float, nargs=2, required=False, default=[-1, 1], help="Range for accelerations (min max)")
    parser.add_argument("-m", "--mass_range", type=float, nargs=2, required=True, help="Range for masses (min max)")

    # Add a flag to set all initial accelerations to zero
    parser.add_argument("--zero_acc", type=bool, required=False, default=True, help="Set accelerations to zero")

    # Add an argument to specify the name of the output CSV file
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file name")

    # Parse the arguments from the command line
    args = parser.parse_args()

    # Call the function to generate initial conditions using the provided arguments
    generate_initial_conditions(
        num_bodies=args.num_bodies,
        pos_range=tuple(args.pos_range),
        vel_range=tuple(args.vel_range),
        acc_range=tuple(args.acc_range),
        mass_range=tuple(args.mass_range),
        zero_acc=args.zero_acc,
        output_file=args.output
    )
