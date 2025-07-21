import pyvista as pv
import argparse

def read_vtk(filename):
    # Load the .vtk file
    mesh = pv.read(filename)

    # Print basic info
    print(f"Loaded: {filename}")
    print(f"Number of points: {mesh.n_points}")
    print(f"Point coordinates:\n{mesh.points[:5]} ...")
    
    # Print available scalar fields
    if mesh.point_data:
        print("Available scalar fields:")
        for name in mesh.point_data:
            print(f" - {name}: min={mesh[name].min()}, max={mesh[name].max()}")

    # Optional: Plot
    mesh.plot(scalars=list(mesh.point_data.keys())[0] if mesh.point_data else None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read and visualize surface .vtk file.")
    parser.add_argument("filename", help="Path to the .vtk file")
    args = parser.parse_args()

    read_vtk(args.filename)

