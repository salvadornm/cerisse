import numpy as np
import matplotlib.pyplot as plt

def generate_circle_dat(
    center=(0.0, 0.0),
    radius=0.1,
    n_points=200,
    filename="circle.dat",
    plot=True
):
    cx, cy = center

    theta = np.linspace(0, 2*np.pi, n_points, endpoint=False)
    x = cx + radius * np.cos(theta)
    y = cy + radius * np.sin(theta)

    data = np.column_stack((x, y))
    np.savetxt(filename, data, fmt="%.10f")

    print(f"Circle data saved to: {filename}")
    print(f"Center = {center}, Radius = {radius}, Points = {n_points}")

    if plot:
        plt.figure()
        plt.plot(x, y, 'o-', markersize=2)
        plt.gca().set_aspect('equal')
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Generated Circle")
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    generate_circle_dat(
        center=(0.7, 0.0),   # Circle Center
        radius=0.1,          # Radius
        n_points=1000,       # Number of Points
        filename="circle7.dat",
        plot=True
    )

