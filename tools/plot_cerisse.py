
import yt
import sys
import os
import glob
import numpy as np
import argparse

# Define the Mach number field
def _mach_number(field, data):
    gamma = 1.4 # Ratio of specific heats
    p = data['boxlib', 'pressure']
    rho = data['boxlib', 'Density']
    # Ensure positive pressure and density to avoid sqrt errors
    p = np.maximum(p, 1e-10)
    rho = np.maximum(rho, 1e-10)
    c = np.sqrt(gamma * p / rho)
    v_mag = data['gas', 'velocity_magnitude']
    return v_mag / c

def plot_field(plotfile, field_name, log_scale=False, cmap='jet'):
    """
    Plots a specific field from a cerisse/AMReX plotfile using yt.
    """
    if not os.path.exists(plotfile):
        print(f"Error: Plotfile {plotfile} does not exist.")
        return

    ds = yt.load(plotfile)

    # Add Mach number field if needed
    if field_name == 'Mach':
        ds.add_field(('gas', 'mach_number'), function=_mach_number, sampling_type='cell', units='')
        target_field = ('gas', 'mach_number')
    elif field_name.lower() == 'pressure':
        # Handle case sensitivity if needed, but usually it's ('boxlib', 'pressure')
        # We'll try to find the exact match in the dataset or default to boxlib
        target_field = ('boxlib', 'pressure')
    elif field_name in ['Density', 'density']:
        target_field = ('boxlib', 'Density')
    else:
        # Assume the user passed a valid field name, possibly needing a tuple prefix
        # If it's a simple string, yt often finds it. If not, we might need to guess 'boxlib'
        target_field = field_name

    print(f"Plotting {field_name}...")

    # Define a masked field generator
    def _masked_field(field, data):
        # Get the original field data
        if field_name == 'Mach':
            # Re-calculate Mach locally
            gamma = 1.4
            p = data['boxlib', 'pressure']
            rho = data['boxlib', 'Density']
            p = np.maximum(p, 1e-10)
            rho = np.maximum(rho, 1e-10)
            c = np.sqrt(gamma * p / rho)
            val = data['gas', 'velocity_magnitude'] / c
        else:
            val = data[target_field]
            
        # Get sld
        sld = data['boxlib', 'sld']
        
        # Create mask: where sld > 0.5, set to NaN
        mask = sld > 0.5
        
        # Return masked array
        val_masked = np.copy(val)
        val_masked[mask] = np.nan
        return val_masked

    # Register the masked field
    # We use a unique name to avoid conflicts
    # If target_field is a tuple, use the second part for the name
    base_name = target_field[1] if isinstance(target_field, tuple) else target_field
    masked_field_name = ('gas', f'masked_{base_name}')
    
    ds.add_field(masked_field_name, function=_masked_field, sampling_type='cell', units='', force_override=True)

    # Create plot with masked field
    try:
        slc = yt.SlicePlot(ds, 'z', masked_field_name)
        slc.set_log(masked_field_name, log_scale)
        slc.set_cmap(masked_field_name, cmap)
        
        # Set NaN color to white
        slc.set_background_color(masked_field_name, 'white')
        
        # Adjust title
        slc.annotate_title(f"Slice of {field_name}")

        # Save
        saved_files = slc.save()
        print(f"Saved plot to: {saved_files}")
    except Exception as e:
        print(f"Error plotting {field_name}: {e}")
        print("Available fields:", ds.field_list)

def main():
    parser = argparse.ArgumentParser(description="Plot fields from Cerisse/AMReX plotfiles.")
    parser.add_argument("plotfile", nargs="?", help="Path to the plotfile")
    parser.add_argument("--plot", action="store_true", help="Enable plotting")
    parser.add_argument("--field", help="Name of the field to plot (e.g., 'pressure', 'Mach', 'Density')")
    
    args = parser.parse_args()

    # Determine plotfile path
    plotfile = args.plotfile
    
    if not plotfile:
        # Determine base_dir relative to this script's location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.join(script_dir, "../wrk/retro/plot")
        
        search_pattern = os.path.join(base_dir, "plt*")
        files = sorted(glob.glob(search_pattern))
        # Filter out .temp files
        files = [f for f in files if not f.endswith('.temp')]
        
        if files:
            plotfile = files[-1]
            print(f"No file specified. Using latest found: {plotfile}")
        else:
            print(f"No plotfiles found in default location. Please specify a file.")
            return

    if args.plot:
        if args.field:
            # Determine log scale based on field name
            # Default to False, set True for pressure
            is_log = False
            if args.field.lower() == 'pressure':
                is_log = True
            
            plot_field(plotfile, args.field, log_scale=is_log, cmap='jet')
        else:
            # Default behavior: Plot Pressure and Mach
            print("No field specified with --field. Plotting default fields (Pressure, Mach).")
            plot_field(plotfile, 'pressure', log_scale=True, cmap='jet')
            plot_field(plotfile, 'Mach', log_scale=False, cmap='jet')
    else:
        # If --plot is not specified, maybe just list fields or do nothing?
        # For now, let's assume if the user runs the script they probably want to plot 
        # even if they forgot --plot, OR we strictly follow the flag.
        # Given the user request "python3 ... --plot --field ...", I will strictly follow --plot
        # but if no arguments are provided at all (except maybe filename), I'll default to plotting defaults for convenience.
        if args.field:
             print("Please add --plot to generate the plot.")
        else:
             # Legacy behavior: if no flags, plot defaults
             print("No flags specified. Plotting default fields (Pressure, Mach).")
             plot_field(plotfile, 'pressure', log_scale=True, cmap='jet')
             plot_field(plotfile, 'Mach', log_scale=False, cmap='jet')

if __name__ == "__main__":
    main()
