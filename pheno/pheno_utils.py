#!/usr/bin/env python
"""Some data rejiggling functions"""

import numpy as np
from osgeo import gdal

def resample_dataset ( fname, x_factor, y_factor, method="mean", \
            data_min=-1000, data_max=10000 ):
    """This function resamples a GDAL dataset (single band) by a factor of
    (``x_factor``, ``y_factor``) in x and y. By default, the only method used
    is to calculate the mean. The ``data_min`` and ``data_max`` parameters are
    used to mask out pixels in value"""

    # First open the file
    fname = 'HDF4_EOS:EOS_GRID:"%s":' % fname + \
            'MOD_Grid_monthly_CMG_VI:CMG 0.05 Deg Monthly NDVI'
    gdal_data = gdal.Open ( fname )
    # Get raster sizes
    nx = gdal_data.RasterXSize
    ny = gdal_data.RasterYSize
    # Calculate output raster size
    nnx = nx/x_factor
    nny = ny/y_factor
    # Reshape the raster data...
    B = np.reshape ( gdal_data.ReadAsArray(), ( nnx, x_factor, nny, y_factor ) )
    B = np.ma.array ( B, mask=np.logical_or ( B <= data_min, B >= data_max) )
    # Re-jiggle the dimensions so we can easily average over then
    C = np.transpose ( B, (0, 2, 1, 3 ) )
    if method == "mean":
        reduced_raster = np.mean ( np.mean ( C, axis=-1), axis=-1 )
    else:
        raise NotImplemented, "Only mean reduction supported by now"
    return reduced_raster

def save_raster ( fname_out, raster_in, cell_size, \
        driver="GTiff",dtype=gdal.GDT_Float32 ):
    """This function saves a raster to a filename. The raster must either be
    two-dimensional, or three dimensional, with the first dimension being the
    number of bands. By default, we use GeoTIFF output."""
    
    drv = gdal.CreateDriverByName ( driver )
    # Get shapes
    try:
        ( n_bands, nx, ny ) = raster_in.shape
    except ValueError:
        ( nx, ny ) = raster_in.shape
        n_bands = 1
    # Create output file
    dst_ds = drv.Create ( fname_out, nx, ny, n_bands, dtype )
    dst_ds.SetGeoTransform( [-0.75, cell_size, 0.0, 90.75, 0.0, -cell_size])
    dst_ds.SetProjection ( 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84"' + \
    ',6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],' + \
    'PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",' + \
    '0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]' )

    for b in xrange ( n_bands ):
        try:
            dst_ds.GetRasterBand ( b+1 ).WriteArray ( \
                raster_in [ b, :, :].astype(np.float32) )
        except IndexError:
            dst_ds.GetRasterBand ( b+1 ).WriteArray ( \
                raster_in [ :, :].astype(np.float32) )
    dst_ds = None

def process_vi_files ( data_dir, fname_out, cell_size=1.5, vi="NDVI" ):
    """This function scans all the MODIS HDF files, and process them in annual
    chunks"""
    import glob
    # This globs all the files ending in HDF. I'm assuming that that's the ones
    # I want. Ok, will also do the M?D13C2 bit too...
    files = glob.glob ( "%s/M*D13C2.*.hdf" % data_dir )
    files.sort()
    files = np.array ( files )
    years = np.array( [int(s.split(".")[1][1:5]) for s in files] )
    nny = 180./cell_size
    nnx = 360./cell_size
    x_factor = cell_size/0.05 # CMG cell size is 0.05 degrees
    y_factor = cell_size/0.05 # CMG cell size is 0.05 degrees
    
    for y in np.unique ( years ):
        print "Doing year ", y
        year_sel = ( years == y )
        annual = np.zeros ( ( 12, nny, nnx ) )
        for ( i, f_in ) in enumerate ( files[ year_sel ] ):
            annual [i, :, : ] = resample_dataset ( f_in, x_factor, y_factor )
        save_raster ( "%s_%04d.tif" % ( fname_out, y ), annual, cell_size )
        print "Saved to %s_%04d.tif" % ( fname_out, y )
            
    print "Finished"