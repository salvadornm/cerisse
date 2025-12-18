import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import sys

def visualize_2d_vtp(filename, field_name=None):
    """
    Visualizes a 2D VTP file (which contains lines/contours).
    Provides two views:
    1. Spatial View: The shape in X-Y plane colored by the scalar field.
    2. Data View: The scalar field value plotted against the point index (or X coordinate).
    """
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)

    # Load the mesh
    mesh = pv.read(filename)
    
    # Check dimensionality
    if mesh.n_points == 0:
        print("Error: Mesh has no points.")
        sys.exit(1)
        
    # Get points (ignore Z for 2D)
    x = mesh.points[:, 0]
    y = mesh.points[:, 1]
    
    # Determine which field to plot
    if field_name is None:
        if mesh.cell_data:
            field_name = list(mesh.cell_data.keys())[0]
        elif mesh.point_data:
            field_name = list(mesh.point_data.keys())[0]
        else:
            print("Error: No scalar data found in file.")
            sys.exit(1)
    
    print(f"Visualizing field: '{field_name}'")
    
    # Extract data
    # Note: Cerisse eib.h writes to CellData (one value per line segment)
    # We might want to interpolate to points for smoother plotting, 
    # or just plot cell centers.
    
    data = None
    is_cell_data = False
    
    if field_name in mesh.cell_data:
        # Convert Cell Data to Point Data for better visualization
        # This avoids issues with cell_centers() on some mesh types
        print("Converting Cell Data to Point Data...")
        mesh = mesh.cell_data_to_point_data()
        data = mesh.point_data[field_name]
        
        # Use point coordinates
        plot_x_coord = mesh.points[:, 0]
        plot_y_coord = mesh.points[:, 1]
        
    elif field_name in mesh.point_data:
        data = mesh.point_data[field_name]
        plot_x_coord = x
        plot_y_coord = y
    else:
        print(f"Error: Field '{field_name}' not found.")
        print(f"Available Cell Data: {list(mesh.cell_data.keys())}")
        print(f"Available Point Data: {list(mesh.point_data.keys())}")
        sys.exit(1)

    # --- Visualization ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Spatial View (Shape colored by value)
    # We use a scatter plot here as it's easiest to handle general point clouds
    sc = axes[0].scatter(plot_x_coord, plot_y_coord, c=data, cmap='jet', s=20)
    axes[0].set_title(f"Spatial Distribution: {field_name}")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    axes[0].axis('equal')
    axes[0].grid(True, linestyle='--', alpha=0.6)
    plt.colorbar(sc, ax=axes[0], label=field_name)

    # Plot 2: Data Distribution
    # We plot against X-coordinate. 
    # (If the shape is vertical, plotting against Y might be better, but X is standard for airfoils etc)
    # Sort by X to make the line plot make sense if it's a function of X
    # If it's a closed loop (like a cylinder), plotting vs Index (Arc Length) is often better.
    
    # Let's try plotting vs Index first (Arc Length proxy)
    indices = np.arange(len(data))
    axes[1].plot(indices, data, 'o-', markersize=4, label='Data')
    axes[1].set_title(f"{field_name} vs Index (Along Surface)")
    axes[1].set_xlabel("Index (approx arc length)")
    axes[1].set_ylabel(field_name)
    axes[1].grid(True, linestyle='--', alpha=0.6)
    
    # Optional: If it looks like a function of X (e.g. top and bottom surfaces), 
    # plotting vs X is useful. Let's add a second line for that if requested, 
    # but for now Index is safer for general shapes.

    plt.tight_layout()
    
    output_img = filename + ".png"
    plt.savefig(output_img)
    print(f"Plot saved to: {output_img}")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize 2D VTP file (Matplotlib).")
    parser.add_argument("filename", help="Path to the .vtp file")
    parser.add_argument("--field", type=str, help="Scalar field to plot", default=None)
    
    args = parser.parse_args()
    
    visualize_2d_vtp(args.filename, args.field)
