#Import Libraries
import pickle 
import matplotlib.pyplot as plt
from pathlib import Path
import os
from datetime import datetime

#Add a data path to the folder with the data in pkl files
data_path = r'C:\Physics REU\BAND\seeback_data_new'

#Creating a folder to hold the figures - choice between saving in current folder or home folder
fig_path = Path.home() / "seebeck.figs"
#fig_path = Path.cwd() / "seebeck.figs"
fig_path.mkdir(parents=True, exist_ok=True)
data_dicts = {}

# Loops through all files in the folder
for filename in os.listdir(data_path):
    if filename.endswith('.pkl'):
        file_path = os.path.join(data_path, filename)
        try:
            with open(file_path, 'rb') as f:
                data = pickle.load(f)
                data_dicts[filename] = data
                print(f"Loaded {filename} successfully.")
        except Exception as e:
            print(f"Error loading {filename}: {e}")

#print(data_dicts.keys())


def generate_colors(n, cmap_name='tab10'):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / n) for i in range(n)]



### Create 2x2 Graph w Seebeck_i and DOS(x,y) by bin size for each file ###

for filename in data_dicts: 
    colors_individual = generate_colors(len(data_dicts[filename].keys()))
    fig, axes = plt.subplots(2, 2, figsize=(15, 5)) 
    for idx, bin_key in enumerate(data_dicts[filename].keys()):
        items = dict(data_dicts[filename][bin_key].items())
        temp_array = items["temp"]
        s_xx = items["seebeck_xx"]
        s_yy = items["seebeck_yy"]
        s_zz = items["seebeck_zz"]
        dos_x = items["dos_x"]
        dos_y = items["dos_y"]

        color = colors_individual[idx]

        axes[0,0].plot(temp_array, s_xx, label = f'Bin{bin_key}', color=color)
        axes[0,1].plot(temp_array, s_yy, label = f'Bin{bin_key}', color=color)
        axes[1,0].plot(temp_array, s_zz, label = f'Bin{bin_key}', color=color)
        axes[1,1].plot(dos_x, dos_y, label = f'Bin{bin_key}', color=color)
        #print(items)
        #print("End of entry ----------------------------------")

    axes[0,0].set_title('Seebeck XX')
    axes[0,1].set_title('Seebeck YY')
    axes[1,0].set_title('Seebeck ZZ')
    axes[1,1].set_title('Density of States')

    for idx, ax in enumerate(axes.flatten()):
        if idx < 3:
             ax.set_xlabel('Temperature (K)')
             ax.set_ylabel('Seebeck (µV/K)')
                
        else:
            ax.set_xlabel('DOS (x)')
            ax.set_ylabel('DOS (y)')
        ax.legend(loc='best')
    
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    file_path = os.path.join(fig_path, f'{filename}_{timestamp}.png')
    plt.savefig(file_path, dpi=300)
    plt.close(fig)


### Create a 2x2 graph w Seebeck_i and DOS(x,y) by mesh size for highest bins
mesh_fig, mesh_axes = plt.subplots(2, 2, figsize=(15, 5)) 
colors = generate_colors(data_dicts.__len__())
largest_bins_dict = {}

for filename, bins in data_dicts.items():
    max_npts = max(bins.keys())  
    largest_bin_data = bins[max_npts]

    nkpt = largest_bin_data['nkpt']
    largest_bins_dict[nkpt] = largest_bin_data

for idx, mesh_key in enumerate(largest_bins_dict):
    items = dict(largest_bins_dict[mesh_key].items())
    temp_array = items["temp"]
    s_xx = items["seebeck_xx"]
    s_yy = items["seebeck_yy"]
    s_zz = items["seebeck_zz"]
    dos_x = items["dos_x"]
    dos_y = items["dos_y"]
    color = colors[idx]

    mesh_axes[0,0].plot(temp_array, s_xx, label = f'Mesh{mesh_key}', color=color)
    mesh_axes[0,1].plot(temp_array, s_yy, label = f'Mesh{mesh_key}', color=color)
    mesh_axes[1,0].plot(temp_array, s_zz, label = f'Mesh{mesh_key}', color=color)
    mesh_axes[1,1].plot(dos_x, dos_y, label = f'Mesh{mesh_key}', color=color)


# Customize and Plot the Mesh Graph
mesh_axes[0,0].set_title('Seebeck XX')
mesh_axes[0,1].set_title('Seebeck YY')
mesh_axes[1,0].set_title('Seebeck ZZ')
mesh_axes[1,1].set_title('Density of States')

for idx, ax in enumerate(mesh_axes.flatten()):
        if idx < 3:
             ax.set_xlabel('Temperature (K)')
             ax.set_ylabel('Seebeck (µV/K)')
                
        else:
            ax.set_xlabel('DOS (x)')
            ax.set_ylabel('DOS (y)')
        ax.legend(loc='best')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
mesh_file_path = os.path.join(fig_path, f'mesh_graph_{timestamp}.png')
plt.savefig(mesh_file_path, dpi=300)
plt.close(mesh_fig)



