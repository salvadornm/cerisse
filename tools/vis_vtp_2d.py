"""
vis_vtp_2d.py

Description:
    This script is designed specifically for visualizing 2D surface data (contours/lines)
    from Cerisse .vtp files. 
    
    It provides two modes:
    1. Info Mode (Default): Prints mesh statistics and available data fields.
    2. Plot Mode (--plot): Generates a static Matplotlib figure.
       - For each file, it creates a row with two subplots:
         Left: Spatial View (Shape colored by value)
         Right: Data View (Value vs Index/Arc Length)

Usage:
    1. Inspect file information only:
       python vis_vtp_2d.py path/to/file.vtp

    2. Visualize one or more files (saves to .png automatically):
       python vis_vtp_2d.py file1.vtp file2.vtp --plot --field Pressure

Dependencies:
    - pyvista
    - matplotlib
    - numpy
"""

import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys
import math

def print_file_info(filename):
    """
    Reads a .vtp file and prints information about it (similar to vis_vtp_3d.py).
    """
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        return None

    try:
        mesh = pv.read(filename)
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None

    print(f"Successfully loaded: {filename}")
    print("-" * 40)
    
    # Mesh statistics
    print(f"Number of Points: {mesh.n_points}")
    print(f"Number of Cells:  {mesh.n_cells}")
    
    # Check for Cell Data
    if mesh.cell_data:
        print("\nAvailable Cell Data (Scalars):")
        for name in mesh.cell_data:
            data = mesh.cell_data[name]
            print(f" - {name:<15}: Min = {data.min():.4e}, Max = {data.max():.4e}")
    else:
        print("\nNo Cell Data found.")

    # Check for Point Data
    if mesh.point_data:
        print("\nAvailable Point Data:")
        for name in mesh.point_data:
            data = mesh.point_data[name]
            print(f" - {name:<15}: Min = {data.min():.4e}, Max = {data.max():.4e}")

    print("-" * 40)
    return mesh

def load_and_extract_data(filename, field_name=None):
    """
    Helper function to load mesh and extract data.
    Returns: (plot_x_coord, plot_y_coord, data, actual_field_name)
    """
    if not os.path.exists(filename):
        return None, None, None, None

    try:
        mesh = pv.read(filename)
    except:
        return None, None, None, None
    
    if mesh.n_points == 0:
        return None, None, None, None

    # Determine field
    if field_name is None:
        if mesh.cell_data:
            field_name = list(mesh.cell_data.keys())[0]
        elif mesh.point_data:
            field_name = list(mesh.point_data.keys())[0]
        else:
            return None, None, None, None

    data = None
    plot_x = None
    plot_y = None

    if field_name in mesh.cell_data:
        # Convert Cell Data to Point Data
        mesh = mesh.cell_data_to_point_data()
        data = mesh.point_data[field_name]
        plot_x = mesh.points[:, 0]
        plot_y = mesh.points[:, 1]
    elif field_name in mesh.point_data:
        data = mesh.point_data[field_name]
        plot_x = mesh.points[:, 0]
        plot_y = mesh.points[:, 1]
    else:
        print(f"Warning: Field '{field_name}' not found in {filename}. Skipping.")
        return None, None, None, None
    
    # Ensure data is a flat numpy array for Matplotlib
    if data is not None:
        data = np.array(data)
        if data.ndim > 1:
            data = data.flatten()
            
        # Handle NaNs: Replace with 0 or min value to avoid plotting errors
        if np.isnan(data).any():
            print(f"Warning: Data for '{field_name}' contains NaNs. Replacing with 0 for visualization.")
            data = np.nan_to_num(data, nan=0.0)
        
    return plot_x, plot_y, data, field_name

def visualize_2d_vtp(filenames, field_names=None):
    if not filenames:
        print("No files provided.")
        return

    # Prepare list of (filename, field) to plot
    plot_tasks = []
    for filename in filenames:
        if not field_names:
            plot_tasks.append((filename, None))
        else:
            for field in field_names:
                plot_tasks.append((filename, field))

    n_plots = len(plot_tasks)
    print(f"Generating {n_plots} plots...")
    
    # Layout: One row per task, 2 columns (Spatial, Data)
    n_rows = n_plots
    n_cols = 2
    
    # Adjust figure size based on number of rows
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 5 * n_rows))
    
    # Ensure axes is always 2D array [row, col]
    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)
        
    for i, (filename, req_field) in enumerate(plot_tasks):
        x, y, data, fname = load_and_extract_data(filename, req_field)
        
        ax_spatial = axes[i, 0]
        ax_data = axes[i, 1]
        
        if x is not None:
            # Plot 1: Spatial
            sc = ax_spatial.scatter(x, y, c=data, cmap='jet', s=20)
            ax_spatial.set_title(f"{os.path.basename(filename)}\nSpatial: {fname}", fontsize=10)
            ax_spatial.set_xlabel("X")
            ax_spatial.set_ylabel("Y")
            ax_spatial.axis('equal')
            ax_spatial.grid(True, linestyle='--', alpha=0.6)
            plt.colorbar(sc, ax=ax_spatial, label=fname)

            # Plot 2: Data
            indices = np.arange(len(data))
            ax_data.plot(indices, data, 'o-', markersize=4, label='Data')
            ax_data.set_title(f"{os.path.basename(filename)}\n{fname} vs Index", fontsize=10)
            ax_data.set_xlabel("Index")
            ax_data.set_ylabel(fname)
            ax_data.grid(True, linestyle='--', alpha=0.6)
        else:
            ax_spatial.text(0.5, 0.5, f"Error loading {os.path.basename(filename)}", ha='center', va='center')
            ax_data.axis('off')

    plt.tight_layout()
    
    # Determine output filename
    if n_plots == 1:
        output_img = filenames[0] + ".png"
    else:
        output_img = "comparison_plot_2d.png"
        
    plt.savefig(output_img)
    print(f"Plot saved to: {output_img}")
    # plt.show() # Optional: Comment out if running on headless server

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize 2D VTP file(s).")
    parser.add_argument("filenames", nargs='+', help="Path to the .vtp file(s)")
    parser.add_argument("--plot", action="store_true", help="Generate visualization plots")
    parser.add_argument("--field", nargs='+', type=str, help="Scalar field(s) to plot", default=None)
    
    args = parser.parse_args()
    
    # Always print info for all files first
    for f in args.filenames:
        print_file_info(f)

    # If plot flag is set, generate plots
    if args.plot:
        visualize_2d_vtp(args.filenames, args.field)

