from RadFiled3D.RadFiled3D import FieldStore, CartesianRadiationField, HistogramVoxel, RadiationFieldMetadataV1
import argparse
import numpy as np
import os
import plotly.graph_objects as go
import plotly
import plotly.subplots


def get_kerma_field_from(field: CartesianRadiationField) -> np.ndarray:
    for channel in field.get_channel_names():
        for layer in field.get_channel(channel).get_layers():
            if layer == 'kerma_air':
                print(f"Found kerma_air in field with unit {field.get_channel(channel).get_layer_unit(layer)}")
                return field.get_channel(channel).get_layer_as_ndarray(layer)   

    if "xray_beam" in field.get_channel_names():
        print("No kerma_air field found in field1, using xray_beam instead")
        kerma_field = field.get_channel("xray_beam").get_layer_as_ndarray("spectrum") * field.get_channel("xray_beam").get_layer_as_ndarray("hits")[:, :, :, np.newaxis]
        hist_voxel: HistogramVoxel = field.get_channel("xray_beam").get_voxel_flat("spectrum", 0)
        bin_width = hist_voxel.get_histogram_bin_width()
        bin_edges = np.linspace(0.0, hist_voxel.get_bins() * bin_width, hist_voxel.get_bins())
        
        # expand bin edges histogram to match the shape of the field [x,y,z,bins]
        bin_edges = bin_edges[np.newaxis, np.newaxis, np.newaxis, :]
        bin_edges = np.repeat(bin_edges, kerma_field.shape[0], axis=0)
        bin_edges = np.repeat(bin_edges, kerma_field.shape[1], axis=1)
        bin_edges = np.repeat(bin_edges, kerma_field.shape[2], axis=2)

        kerma_field *= bin_edges
        kerma_field = kerma_field.sum(axis=-1)
    else:
        raise ValueError("No kerma_air field found in field")
    return kerma_field


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two radiation fields')
    parser.add_argument('--file1', type=str, help='First radiation field file')
    parser.add_argument('--file2', type=str, help='Second radiation field file')
    parser.add_argument('--flip_one_around_z', action='store_true', help='Flip the first field around the z-axis')

    args = parser.parse_args()
    
    should_flip_around_z = args.flip_one_around_z
    field1: CartesianRadiationField = FieldStore.load(args.file1)
    field2: CartesianRadiationField = FieldStore.load(args.file2)

    metadata1: RadiationFieldMetadataV1 = FieldStore.load_metadata(args.file1)
    metadata2: RadiationFieldMetadataV1 = FieldStore.load_metadata(args.file2)

    fig = go.Figure()
    # plot the source spectra
    if "tube_spectrum" in metadata1.get_dynamic_metadata_keys():
        spectrum: HistogramVoxel = metadata1.get_dynamic_metadata("tube_spectrum")
        spectrum.normalize()
        bin_width = spectrum.get_histogram_bin_width() / 1000.0
        bin_edges = np.linspace(0.0, spectrum.get_bins() * bin_width, spectrum.get_bins())
        fig.add_trace(go.Scatter(x=bin_edges, y=spectrum.get_histogram(), mode='lines', name="tube spectrum: " + os.path.basename(args.file1)))
    if "tube_spectrum" in metadata2.get_dynamic_metadata_keys():
        spectrum: HistogramVoxel = metadata2.get_dynamic_metadata("tube_spectrum")
        spectrum.normalize()
        bin_width = spectrum.get_histogram_bin_width() / 1000.0
        bin_edges = np.linspace(0.0, spectrum.get_bins() * bin_width, spectrum.get_bins())
        fig.add_trace(go.Scatter(x=bin_edges, y=spectrum.get_histogram(), mode='lines', name="tube spectrum: " + os.path.basename(args.file2)))
    fig.update_xaxes(title_text="Energy (keV)")
    fig.update_yaxes(title_text="Probability density")
    fig.update_layout(title_text="Source spectra")
    fig.show()

    field_dim = field1.get_field_dimensions()
    assert field_dim == field2.get_field_dimensions(), f"Field dimensions do not match {field_dim} != {field2.get_field_dimensions()}"
    print(f"Field dimensions match: {field_dim}")

    kerma_field1 = get_kerma_field_from(field1)
    kerma_field2 = get_kerma_field_from(field2)

    kerma_field1 /= kerma_field1.max()
    kerma_field2 /= kerma_field2.max()

    if should_flip_around_z:
        kerma_field1 = np.flip(kerma_field1, axis=0)

    # Side view
    #kerma_field1 = kerma_field1.transpose((0, 2, 1))
    #kerma_field2 = kerma_field2.transpose((0, 2, 1))
    #kerma_field1 = kerma_field1.transpose((2, 0, 1))
    #kerma_field2 = kerma_field2.transpose((2, 0, 1))

    difference = kerma_field1 - kerma_field2
    print(f"Mean difference: {np.mean(difference)}")
    
    relative_derivation_from_field1 = difference[kerma_field1 != 0.0] / kerma_field1[kerma_field1 != 0.0] * 100
    print(f"Mean relative derivation from field1: {np.mean(relative_derivation_from_field1)}%")
    print(f"Median relative derivation from field1: {np.median(relative_derivation_from_field1)}%")
    

    x_steps = np.linspace(0, field_dim.x, kerma_field1.shape[0])
    y_steps = np.linspace(0, field_dim.y, kerma_field1.shape[1])
    z_steps = np.linspace(0, field_dim.z, kerma_field1.shape[2])


    # make subplots and plot field along the planes: xy, xz, yz for field1 only. Plot as 2d images
    for kerma_field, file_name in [(kerma_field1, os.path.basename(args.file1)), (kerma_field2, os.path.basename(args.file2))]:
        fig = plotly.subplots.make_subplots(rows=1, cols=3, subplot_titles=("xy", "xz", "yz"))
        fig.add_trace(go.Heatmap(z=kerma_field[:, :, kerma_field.shape[2]//2], colorscale='Viridis', showscale=False, x=x_steps, y=y_steps), row=1, col=1)
        fig.add_trace(go.Heatmap(z=kerma_field[:, kerma_field.shape[1]//2, :], colorscale='Viridis', showscale=False, x=x_steps, y=z_steps), row=1, col=2)
        fig.add_trace(go.Heatmap(z=kerma_field[kerma_field.shape[0]//2, :, :], colorscale='Viridis', showscale=False, x=y_steps, y=z_steps), row=1, col=3)
        
        fig.update_xaxes(title_text="X (meters)", row=1, col=1)
        fig.update_yaxes(title_text="Y (meters)", row=1, col=1)
        fig.update_xaxes(title_text="X (meters)", row=1, col=2)
        fig.update_yaxes(title_text="Z (meters)", row=1, col=2)
        fig.update_xaxes(title_text="Y (meters)", row=1, col=3)
        fig.update_yaxes(title_text="Z (meters)", row=1, col=3)
        fig.update_layout(title_text=f"{file_name}")
        fig.show()

    # plot both field as 3dsurfaces in the same figure
    fig = go.Figure()
    fig.add_trace(go.Surface(z=kerma_field1[:, kerma_field1.shape[1]//2, :], colorscale='Viridis', showscale=False, x=x_steps, y=y_steps, name=os.path.basename(args.file1)))
    fig.add_trace(go.Surface(z=kerma_field2[:, kerma_field2.shape[1]//2, :], colorscale='aggrnyl', showscale=False, x=x_steps, y=y_steps, name=os.path.basename(args.file2)))
    fig.update_xaxes(title_text="X (meters)")
    fig.update_yaxes(title_text="Y (meters)")
    fig.update_layout(title_text="Field comparison")
    fig.show()
