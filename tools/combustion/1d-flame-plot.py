import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_flame_structure(filename):
    # Load data
    df = pd.read_csv(filename)

    # Extract relevant data
    x = df["X"]
    temperature = df["temp"]
    x_oh = df["X_OH"]

    # Compute temperature gradient
    temp_gradient = np.gradient(temperature, x)
    # max temperature gradient
    max_gradT = np.max(np.abs(temp_gradient))

    flame_front_index = np.argmax(np.abs(temp_gradient))
    flame_front_x = x[flame_front_index]

    # Shift x so flame front is at 0
    shifted_x = x - flame_front_x

    # Estimate flame thickness
    # Define thickness as distance between 10% and 90% of max temperature rise
    T_min = np.min(temperature)
    T_max = np.max(temperature)
    T_10 = T_min + 0.1 * (T_max - T_min)
    T_90 = T_min + 0.9 * (T_max - T_min)

    idx_10 = np.argmin(np.abs(temperature - T_10))
    idx_90 = np.argmin(np.abs(temperature - T_90))
    flame_thickness = np.abs(x[idx_90] - x[idx_10])

    flame_thickness_maxgrad = (T_max - T_min)/max_gradT

    print(f"Flame front position: {flame_front_x:.6f} m")
    print(f"Maximum temperature gradient: {max_gradT:.3e} [K/m]")
    print(f"Estimated Flame Thickness (based on 10-90% temp rise) : {flame_thickness*1000:.6f} [mm]")
    print(f"Estimated Flame Thickness (based on max gradient) : {flame_thickness_maxgrad*1000:.6f} [mm]")

    plotdomain = 0.03 # <----------------------------  adjust for display

    # Plot
    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Temperature profile
    ax1.plot(shifted_x, temperature, color='tab:red', label='Temperature (K)')
    ax1.set_xlabel('Position X (m) (Flame Front at X=0)')
    ax1.set_ylabel('Temperature (K)', color='tab:red')
    ax1.tick_params(axis='y', labelcolor='tab:red')
    ax1.set_xlim(-0.5*plotdomain, 0.5*plotdomain)
    ax1.axvline(0, color='gray', linestyle='--', label='Flame Front')

    # OH mole fraction
    ax2 = ax1.twinx()
    ax2.plot(shifted_x, x_oh, color='tab:blue', label='X_OH')
    ax2.set_ylabel('Mole Fraction of OH', color='tab:blue')
    ax2.tick_params(axis='y', labelcolor='tab:blue')

    # Add flame thickness markers
    #ax1.axvline(x[idx_10] - flame_front_x, color='green', linestyle=':', label='10% Temp Rise')
    #ax1.axvline(x[idx_90] - flame_front_x, color='purple', linestyle=':', label='90% Temp Rise')

    # Legends
    fig.legend(loc="upper right", bbox_to_anchor=(0.9, 0.9))
    plt.title('Flame Structure and Flame Thickness')
    plt.grid(True)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    plot_flame_structure("onedim.csv")

