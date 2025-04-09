import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling

# Open the input raster
with rasterio.open("../R-data/places_fmv_pnas_dryad/1 estimates/places_fmv_vacant.tif") as src:
    # Define target CRS (using EPSG code, e.g., EPSG:4326 for WGS84)
    dst_crs = 'EPSG:2163'
    
    # Calculate the transformation parameters
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)
    
    # Prepare the metadata for the output file
    kwargs = src.meta.copy()
    kwargs.update({
        'crs': dst_crs,
        'transform': transform,
        'width': width,
        'height': height
    })
    
    # Create the output file with new CRS
    with rasterio.open('../R-data/places_fmv_pnas_dryad/1 estimates/output.tif', 'w', **kwargs) as dst:
        # Reproject and save each band
        for i in range(1, src.count + 1):
            reproject(
                source=rasterio.band(src, i),
                destination=rasterio.band(dst, i),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest
            )