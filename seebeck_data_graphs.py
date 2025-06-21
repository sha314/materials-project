#Import Libraries
import pickle 
import matplotlib.pyplot as plt
import os

#Add a path to the folder with the data in pkl files
folder_path = r'C:\Physics REU\BAND'
data_dicts = {}

# Loops through all files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.pkl'):
        file_path = os.path.join(folder_path, filename)
        try:
            with open(file_path, 'rb') as f:
                data = pickle.load(f)
                data_dicts[filename] = data
                print(f"Loaded {filename} successfully.")
        except Exception as e:
            print(f"Error loading {filename}: {e}")

#print(data_dicts)

#Create Subplot Graph
fig, axes = plt.subplots(1, 3, figsize=(15, 5)) 


def generate_colors(n, cmap_name='tab10'):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / n) for i in range(n)]

#Checking if only one file in folder, which assumes graph of varying bin sizes
#If multiple files are present, will plot the highest bin size based on mesh sizes

if data_dicts.__len__() == 1:
    colors = generate_colors(len(data.keys()))
    # Loop Across Bin Sizes to Create Graphs
    for idx, i in enumerate(data.keys()):
        items = dict(data[i].items())
        temp_array = items["temp"]
        s_xx = items["seebeck_xx"]
        s_yy = items["seebeck_yy"]
        s_zz = items["seebeck_zz"]
        color = colors[idx]

        axes[0].plot(temp_array, s_xx, label = f'Bin{i}', color=color)
        axes[1].plot(temp_array, s_yy, label = f'Bin{i}', color=color)
        axes[2].plot(temp_array, s_zz, label = f'Bin{i}', color=color)
elif data_dicts.__len__() > 1:
    colors = generate_colors(data_dicts.__len__())
    largest_bins_dict = {}

    for filename, bins in data_dicts.items():
        max_npts = max(bins.keys())  
        largest_bin_data = bins[max_npts]
    
        nkpt = largest_bin_data['nkpt']
        largest_bins_dict[nkpt] = largest_bin_data
    
    for idx, i in enumerate(largest_bins_dict):
        items = dict(largest_bins_dict[i].items())
        temp_array = items["temp"]
        s_xx = items["seebeck_xx"]
        s_yy = items["seebeck_yy"]
        s_zz = items["seebeck_zz"]
        color = colors[idx]

        axes[0].plot(temp_array, s_xx, label = f'Mesh{i}', color=color)
        axes[1].plot(temp_array, s_yy, label = f'Mesh{i}', color=color)
        axes[2].plot(temp_array, s_zz, label = f'Mesh{i}', color=color)
else:
    print("No files read. Make sure the correct folder path is inputted")


# Customize and Plot the Graph
axes[0].set_title('Seebeck XX')
axes[1].set_title('Seebeck YY')
axes[2].set_title('Seebeck ZZ')

for ax in axes:
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Seebeck (µV/K)')
    ax.legend(loc='best')

plt.tight_layout()
plt.show()


  
