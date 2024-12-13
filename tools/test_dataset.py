from pyRadiationField.pyRadiationField import CartesianRadiationField, FieldStore
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import os


data_dir = os.path.join(os.path.dirname(__file__), "Datasets/120kV-Alderson-C100-Sequence")

for file_name in [f for f in os.listdir(data_dir) if f.endswith(".rf")]:
    file_name = os.path.join(data_dir, file_name)

    print(f"Loading field: {file_name}")

    try:
        field: CartesianRadiationField = FieldStore.load(file_name)
    except Exception as e:
        print(f"Error loading field: {file_name}")
        print(e)
        continue
    meta_data = FieldStore.load_metadata(file_name)
    print(f"Field loaded: {field}")
    print(f"Primay Particles: {meta_data.get_header().simulation.primary_particle_count:e}")

    scatter_component = field.get_channel("scatter_field")
    beam_component = field.get_channel("xray_beam")

    scatter_energy_field = scatter_component.get_layer_as_ndarray("energy")
    spectrum = scatter_component.get_layer_as_ndarray("spectrum")
    spectrum_mean_sum = spectrum.sum(axis=-1)   

    print("scatter field")
    print("Min energy: ", scatter_energy_field.min())
    print("Max energy: ", scatter_energy_field.max())
    print("Mean energy: ", scatter_energy_field[scatter_energy_field != 0.0].mean())
    print("Std energy: ", scatter_energy_field[scatter_energy_field != 0.0].std())
    print("Sum energy: ", scatter_energy_field[scatter_energy_field != 0.0].sum())
    print(f"Relative statistical error scatter: {scatter_component.get_statistical_error('spectrum') * 100}%")
    print(f"Relative statistical error beam   : {beam_component.get_statistical_error('spectrum') * 100}%")

    scatter_error_field = scatter_component.get_layer_as_ndarray("error")
    print("Min error: ", scatter_error_field.min())
    print("Max error: ", scatter_error_field.max())
    print("Mean error: ", scatter_error_field[scatter_error_field != 0.0].mean())
    print("Std error: ", scatter_error_field[scatter_error_field != 0.0].std())
    print("Median error: ", np.median(scatter_error_field[scatter_error_field != 0.0]))

    # plot error field as histogram
    sorted_errors = np.sort(scatter_error_field[scatter_error_field != 0.0].flatten())
    print("Error field 90th percentile: ", sorted_errors[int(0.9 * sorted_errors.size)])
    plt.hist(scatter_error_field[scatter_error_field != 0.0].flatten(), bins=100)
    plt.title("Error field histogram")
    plt.show()

    # plot error field as 3d surface at y=height/2
    x_data, z_data = np.meshgrid(np.arange(scatter_error_field.shape[0]), np.arange(scatter_error_field.shape[2]))
    y_data = scatter_error_field[:, scatter_error_field.shape[1] // 2, :]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("Error field at height/2")
    plt.show()

    print("NaN elements count in scatter field: ", scatter_energy_field[np.isnan(scatter_energy_field)].size)
    print("Hits count in scatter field: ", scatter_component.get_layer_as_ndarray("hits").sum() * meta_data.get_header().simulation.primary_particle_count)

    # calculate the duration of the simulation in seconds
    duration = meta_data.get_dynamic_metadata("simulation_duration_s").get_data()
    print(f"Simulation duration: {duration} seconds")
    # calculate speed in particles per second
    speed = meta_data.get_header().simulation.primary_particle_count / duration
    print(f"Speed: {speed} particles / second")

    z_data = np.sum(scatter_energy_field, axis=2)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(scatter_energy_field.shape[0]), np.arange(scatter_energy_field.shape[1]))

    fig = plt.figure()
    ax = fig.add_subplot(131, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("Scatter field Z-Summed")

    z_data = np.sum(scatter_energy_field, axis=1)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(scatter_energy_field.shape[0]), np.arange(scatter_energy_field.shape[1]))
    ax = fig.add_subplot(132, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("Scatter field Y-Summed")

    z_data = np.sum(scatter_energy_field, axis=0)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(scatter_energy_field.shape[0]), np.arange(scatter_energy_field.shape[1]))
    ax = fig.add_subplot(133, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("Scatter field X-Summed")

    xray_beam_energy_field = beam_component.get_layer_as_ndarray("energy")

    plt.show(block=True)
    fig = plt.figure()

    print("xray beam")
    print("Min energy: ", xray_beam_energy_field.min())
    print("Max energy: ", xray_beam_energy_field.max())
    print("Mean energy: ", xray_beam_energy_field[xray_beam_energy_field != 0.0].mean())
    print("Std energy: ", xray_beam_energy_field[xray_beam_energy_field != 0.0].std())
    print("Sum energy: ", xray_beam_energy_field[xray_beam_energy_field != 0.0].sum())

    print("NaN elements count in xray beam field: ", xray_beam_energy_field[np.isnan(xray_beam_energy_field)].size)
    print("Hits count in xray beam field: ", beam_component.get_layer_as_ndarray("hits").sum() * meta_data.get_header().simulation.primary_particle_count)


    z_data = np.sum(xray_beam_energy_field, axis=2)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(xray_beam_energy_field.shape[0]), np.arange(xray_beam_energy_field.shape[1]))
    ax = fig.add_subplot(131, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("XRay beam field Z-Summed")

    z_data = np.sum(xray_beam_energy_field, axis=1)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(xray_beam_energy_field.shape[0]), np.arange(xray_beam_energy_field.shape[1]))
    ax = fig.add_subplot(132, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("XRay beam field Y-Summed")

    z_data = np.sum(xray_beam_energy_field, axis=0)
    z_data = z_data / z_data.max()
    x_data, y_data = np.meshgrid(np.arange(xray_beam_energy_field.shape[0]), np.arange(xray_beam_energy_field.shape[1]))
    ax = fig.add_subplot(133, projection='3d')
    ax.plot_surface(x_data, y_data, z_data)
    ax.set_title("XRay beam field X-Summed")


    plt.show(block=True)
