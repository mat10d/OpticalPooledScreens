# Features adapted from:
# Bray et al. 2016 Nature Protocols 11:1757-1774
# Cell Painting, a high-content image-based assay for morphological profiling using multiplexed fluorescent dyes
#
# Also helpful: 
# 	https://github.com/carpenterlab/2016_bray_natprot/wiki/What-do-Cell-Painting-features-mean%3F
#	https://raw.githubusercontent.com/wiki/carpenterlab/2016_bray_natprot/attachments/feature_names.txt
#	http://cellprofiler-manual.s3.amazonaws.com/CellProfiler-3.1.5/modules/measurement.html

# Note: their protocol uses 20X objective, 2x2 binning. 
# Length scales needed for feature extraction should technically be correspondingly scaled,
# e.g., 20X with 1x1 binning images should use suggested linear length scales * 2

import numpy as np
from scipy.stats import median_abs_deviation, rankdata # new in version 1.3.0
from scipy.spatial.distance import pdist
from scipy.ndimage.morphology import distance_transform_edt as distance_transform
from scipy.ndimage import map_coordinates
from mahotas.features import haralick, pftas, zernike_moments
from mahotas.thresholding import otsu
from scipy.spatial import ConvexHull
from functools import partial
from itertools import starmap, combinations
from warnings import catch_warnings,simplefilter
import skimage.measure
import skimage.morphology
from skimage import img_as_ubyte
from ops.features import correlate_channels_masked, masked, correlate_channels_all_multichannel
from ops.utils import subimage
from decorator import decorator

## SEGMENTATION

def identify_secondary_objects(image, primary_segmentation, method='propagation', regularization_factor=0.05,
	threshold='otsu', remove_boundary_objects=True):
	if method != 'propagation':
		raise ValueError(f'method={method} not implemented')

	from centrosome.propagate import propagate

	if isinstance(threshold,np.ndarray):
		# pre-thresholded
		thresholded = threshold
	elif isinstance(threshold,(int,float)):
		thresholded = image > threshold
	elif threshold == 'otsu':
		thresholded = image > otsu(image)
	else:
		raise ValueError(f'threshold={threshold} not implemented')
		
	secondary_segmentation,_ = propagate(img_as_ubyte(image),primary_segmentation,thresholded,regularization_factor)

	if remove_boundary_objects:
		cut = np.concatenate([secondary_segmentation[0,:], secondary_segmentation[-1,:],
			secondary_segmentation[:,0], secondary_segmentation[:,-1]])
		secondary_segmentation.flat[np.in1d(secondary_segmentation, np.unique(cut))] = 0

	return secondary_segmentation.astype(np.uint16)

## FEATURES

def apply_extract_features_cp(well_tile,filepattern):
	from ops.io_hdf import read_hdf_image
	from ops.io import read_stack as read
	from ops.firesnake import Snake
	from ops.filenames import name_file as name
	wildcards = {'well':well_tile[0],'tile':well_tile[1]}
	filepattern.update(wildcards)
	stacked = read_hdf_image(name(filepattern))
	nuclei = read(name(filepattern,subdir='process_ph',tag='nuclei',ext='tif'))
	cells = read(name(filepattern,subdir='process_ph',tag='cells',ext='tif'))
	df_result = Snake._extract_phenotype_cp(data_phenotype=stacked,
                                            nuclei=nuclei,
                                            cells=cells,
                                            wildcards=wildcards,
                                            nucleus_channels=[0,1,2,3],
                                            cell_channels=[0,1,2,3],
                                            channel_names=['dapi','tubulin','gh2ax','phalloidin']
                                           )
	df_result.to_csv(name(filepattern,subdir='process_ph',tag='cp_phenotype',ext='csv'))

EDGE_CONNECTIVITY = 2 # cellprofiler uses edge connectivity of 1, which exlucdes pixels catty-corner to a boundary

ZERNIKE_DEGREE = 9

# GRANULARITY_BACKGROUND = 10 #this should be a bit larger than the radius of the features, i.e., "granules", of interest after downsampling
# GRANULARITY_BACKGROUND_DOWNSAMPLE = 1
# GRANULARITY_DOWNSAMPLE = 1
# GRANULARITY_LENGTH = 16

# # MeasureCorrelation:'Measure the intensity correlation between all channels, within all objects.'
# #     Hidden:5
# #     Hidden:3
# #     Select an image to measure:DNA
# #     Select an image to measure:Mito
# #     Select an image to measure:ER
# #     Select an image to measure:RNA
# #     Select an image to measure:AGP
# #     Select where to measure correlation:Within objects
# #     Select an object to measure:Nuclei
# #     Select an object to measure:Cells
# #     Select an object to measure:Cytoplasm

# This module is now called "MeasureColocalization" in CellProfiler

correlation_features ={
	'correlation' : lambda r: [correlate_channels_masked(r,first,second) 
	for first,second in combinations(list(range(r.intensity_image_full.shape[-3])),2)],
	'lstsq_slope' : lambda r: [lstsq_slope(r,first,second) 
	for first,second in combinations(list(range(r.intensity_image_full.shape[-3])),2)],
	# costes threshold algorithm not working well, using otsu threhold instead
	'colocalization' : lambda r: cp_colocalization_all_channels(r,mode='old',threshold='otsu')
}

correlation_features_ch ={
	'correlation' : lambda r, ch1, ch2: correlate_channels_masked(r,ch1,ch2),
	'lstsq_slope' : lambda r, ch1, ch2: lstsq_slope(r,ch1,ch2),
	# costes threshold algorithm not working well, using otsu threshold instead
	'colocalization' : lambda r, ch1, ch2: cp_colocalization(r,ch1,ch2,mode='old',threshold='otsu')
}

correlation_features_multichannel ={
	'correlation' : lambda r: catch_runtime(correlate_channels_all_multichannel)(r),
	'lstsq_slope' : lambda r: lstsq_slope_all_multichannel(r),
	# costes threshold algorithm not working well, using otsu threshold instead
	'colocalization' : lambda r: cp_colocalization_all_channels(r,mode='multichannel',threshold='otsu')
}

correlation_columns = [
'correlation_{first}_{second}',
'lstsq_slope_{first}_{second}'
]

colocalization_columns = [
'overlap_{first}_{second}',
'K_{first}_{second}',
'K_{second}_{first}',
'manders_{first}_{second}',
'manders_{second}_{first}',
'rwc_{first}_{second}',
'rwc_{second}_{first}'
]

colocalization_columns_ch = {
'colocalization_{first}_{second}_0':'overlap_{first}_{second}',
'colocalization_{first}_{second}_1':'K_{first}_{second}',
'colocalization_{first}_{second}_2':'K_{second}_{first}',
'colocalization_{first}_{second}_3':'manders_{first}_{second}',
'colocalization_{first}_{second}_4':'manders_{second}_{first}',
'colocalization_{first}_{second}_5':'rwc_{first}_{second}',
'colocalization_{first}_{second}_6':'rwc_{second}_{first}'
}

correlation_columns_multichannel = {
	'correlation':['correlation_{first}_{second}'],
	'lstsq_slope':['lstsq_slope_{first}_{second}'],
	'colocalization':['overlap_{first}_{second}', 'K_{first}_{second}',
					  'K_{second}_{first}', 'manders_{first}_{second}',
					  'manders_{second}_{first}', 'rwc_{first}_{second}',
					  'rwc_{second}_{first}']
}

# # MeasureGranularity:'Measure the granularity characteristics across all images. Note that this is a per-image measure, and therefore, will not appear in the per-well profiles generated by the per-cell features.'
# #     Image count:5
# #     Object count:1
# #     Select an image to measure:DNA
# #     Subsampling factor for granularity measurements:0.25
# #     Subsampling factor for background reduction:0.25
# #     Radius of structuring element:10
# #     Range of the granular spectrum:16
# #     Select objects to measure:Nuclei
# #     Object count:3
# #     Select an image to measure:Mito
# #     Subsampling factor for granularity measurements:0.25
# #     Subsampling factor for background reduction:0.25
# #     Radius of structuring element:10
# #     Range of the granular spectrum:16
# #     Select objects to measure:Cells
# #     Select objects to measure:Nuclei
# #     Select objects to measure:Cytoplasm
# #     Object count:3
# #     Select an image to measure:ER
# #     Subsampling factor for granularity measurements:0.25
# #     Subsampling factor for background reduction:0.25
# #     Radius of structuring element:10
# #     Range of the granular spectrum:16
# #     Select objects to measure:Cells
# #     Select objects to measure:Nuclei
# #     Select objects to measure:Cytoplasm
# #     Object count:3
# #     Select an image to measure:AGP
# #     Subsampling factor for granularity measurements:0.25
# #     Subsampling factor for background reduction:0.25
# #     Radius of structuring element:10
# #     Range of the granular spectrum:16
# #     Select objects to measure:Cells
# #     Select objects to measure:Cytoplasm
# #     Select objects to measure:Nuclei
# #     Object count:3
# #     Select an image to measure:RNA
# #     Subsampling factor for granularity measurements:0.25
# #     Subsampling factor for background reduction:0.25
# #     Radius of structuring element:10
# #     Range of the granular spectrum:16
# #     Select objects to measure:Cells
# #     Select objects to measure:Cytoplasm
# #     Select objects to measure:Nuclei

# In CellProfiler this is a per-image metric, but is implemented here as a per-object metric.
# to re-produce values from paper, use start_radius = 10, spectrum_length = 16, sample=sample_background=0.25
# values here to optimize for fine speckles in single cells: THESE PARAMETERS ARE HIGHLY EXPERIMENT-DEPENDENT

# granularity_features = {
# 	'granularity_spectrum' : lambda r: granularity_spectrum(r.intensity_image_full, r.image, 
# 		background_radius=GRANULARITY_BACKGROUND, spectrum_length=GRANULARITY_LENGTH, 
# 		downsample=GRANULARITY_DOWNSAMPLE, background_downsample=GRANULARITY_BACKGROUND_DOWNSAMPLE)
# }

# # MeasureObjectIntensity:'Measure the intensity characteristics from all channels, within all objects.'
# # 	Hidden:5
# #     Select an image to measure:DNA
# #     Select an image to measure:ER
# #     Select an image to measure:RNA
# #     Select an image to measure:AGP
# #     Select an image to measure:Mito
# #     Select objects to measure:Nuclei
# #     Select objects to measure:Cytoplasm
# #     Select objects to measure:Cells

intensity_features = {
	'int': lambda r: r.intensity_image[r.image].sum(),
	'mean': lambda r: r.intensity_image[r.image].mean(),
	'std': lambda r: np.std(r.intensity_image[r.image]),
	'max': lambda r: r.intensity_image[r.image].max(),
	'min': lambda r: r.intensity_image[r.image].min(),
	'edge_intensity_feature': lambda r: edge_intensity_features(r.intensity_image,r.filled_image,mode='inner',connectivity=EDGE_CONNECTIVITY),
	'mass_displacement': lambda r: np.sqrt(((np.array(r.local_centroid) - np.array(catch_runtime(lambda r: r.weighted_local_centroid)(r)))**2).sum()),
	'lower_quartile': lambda r: np.percentile(r.intensity_image[r.image],25),
    'median': lambda r: np.median(r.intensity_image[r.image]),
    'mad': lambda r: median_abs_deviation(r.intensity_image[r.image],scale=1),
    'upper_quartile': lambda r: np.percentile(r.intensity_image[r.image],75),
    'center_mass': lambda r: catch_runtime(lambda r: r.weighted_local_centroid)(r), # this property is not cached
    'max_location':lambda r: np.unravel_index(np.argmax(r.intensity_image), (r.image).shape)
    }

intensity_features_ch = {
	'int': lambda r, ch : r.intensity_image_full[ch,r.image].sum(),
	'mean': lambda r, ch: r.intensity_image_full[ch,r.image].mean(),
	'std': lambda r, ch: np.std(r.intensity_image_full[ch,r.image]),
	'max': lambda r, ch: r.intensity_image_full[ch,r.image].max(),
	'min': lambda r, ch: r.intensity_image_full[ch,r.image].min(),
	'edge_intensity_feature': lambda r, ch: edge_intensity_features(r.intensity_image_full[ch],r.filled_image,mode='inner',connectivity=EDGE_CONNECTIVITY),
	'mass_displacement': lambda r, ch: mass_displacement_grayscale(r.local_centroid,r.intensity_image_full[ch]*r.image),
	'lower_quartile': lambda r, ch: np.percentile(r.intensity_image_full[ch,r.image],25),
    'median': lambda r, ch: np.median(r.intensity_image_full[ch,r.image]),
    'mad': lambda r, ch: median_abs_deviation(r.intensity_image_full[ch,r.image],scale=1),
    'upper_quartile': lambda r, ch: np.percentile(r.intensity_image_full[ch,r.image],75),
    'center_mass': lambda r, ch: weighted_local_centroid_grayscale(r.intensity_image_full[ch]*r.image),
    'max_location':lambda r, ch: np.unravel_index(np.argmax(r.intensity_image_full[ch]*r.image), (r.image).shape)
    }

intensity_features_multichannel = {
	'int': lambda r: r.intensity_image[r.image,...].sum(axis=0),
	'mean': lambda r: r.intensity_image[r.image,...].mean(axis=0),
	'std': lambda r: np.std(r.intensity_image[r.image,...],axis=0),
	'max': lambda r: r.intensity_image[r.image,...].max(axis=0),
	'min': lambda r: r.intensity_image[r.image,...].min(axis=0),
	'edge_intensity_feature': lambda r: edge_intensity_features(r.intensity_image,r.filled_image,mode='inner',connectivity=EDGE_CONNECTIVITY),
	'mass_displacement': lambda r: np.sqrt(((np.array(r.local_centroid)[:,None] - catch_runtime(lambda r: r.weighted_local_centroid)(r))**2).sum(axis=0)),
	'lower_quartile': lambda r: np.percentile(r.intensity_image[r.image,...],25,axis=0),
    'median': lambda r: np.median(r.intensity_image[r.image,...],axis=0),
    'mad': lambda r: median_abs_deviation(r.intensity_image[r.image,...],scale=1,axis=0),
    'upper_quartile': lambda r: np.percentile(r.intensity_image[r.image,...],75,axis=0),
    'center_mass': lambda r: catch_runtime(lambda r: r.weighted_local_centroid)(r).flatten(), # this property is not cached
    'max_location':lambda r: np.array(np.unravel_index(np.argmax(r.intensity_image.reshape(-1,*r.intensity_image.shape[2:]),axis=0), (r.image).shape)).flatten()
    }

intensity_columns = {
	'edge_intensity_feature_0':'int_edge',
	'edge_intensity_feature_1':'mean_edge',
	'edge_intensity_feature_2':'std_edge',
	'edge_intensity_feature_3':'max_edge',
	'edge_intensity_feature_4':'min_edge',
	'center_mass_0':'center_mass_r',
	'center_mass_1':'center_mass_c',
	'max_location_0':'max_location_r',
	'max_location_1':'max_location_c'
}

intensity_columns_ch = {
	'{channel}_edge_intensity_feature_0':'{channel}_int_edge',
	'{channel}_edge_intensity_feature_1':'{channel}_mean_edge',
	'{channel}_edge_intensity_feature_2':'{channel}_std_edge',
	'{channel}_edge_intensity_feature_3':'{channel}_max_edge',
	'{channel}_edge_intensity_feature_4':'{channel}_min_edge',
	'{channel}_center_mass_0':'{channel}_center_mass_r',
	'{channel}_center_mass_1':'{channel}_center_mass_c',
	'{channel}_max_location_0':'{channel}_max_location_r',
	'{channel}_max_location_1':'{channel}_max_location_c'
}

intensity_columns_multichannel = {
	'int':['int'],
	'mean':['mean'],
	'std':['std'],
	'max':['max'],
	'min':['min'],
	'mass_displacement':['mass_displacement'],
	'lower_quartile':['lower_quartile'],
	'median':['median'],
	'mad':['mad'],
	'upper_quartile':['upper_quartile'],
	'edge_intensity_feature':['int_edge', 'mean_edge','std_edge', 'max_edge','min_edge'],
	'center_mass':['center_mass_r','center_mass_c'],
	'max_location':['max_location_r','max_location_c']
}

# # MeasureObjectNeighbors:'Measure the adjacency statistics for the cells. Cells within 5 pixels of each other are considered neighbors.'
# #     Select objects to measure:Cells
# #     Select neighboring objects to measure:Cells
# #     Method to determine neighbors:Within a specified distance
# #     Neighbor distance:5
# #     Retain the image of objects colored by numbers of neighbors?:No
# #     Name the output image:ObjectNeighborCount
# #     Select colormap:Default
# #     Retain the image of objects colored by percent of touching pixels?:No
# #     Name the output image:PercentTouching
# #     Select a colormap:Default

# # MeasureObjectNeighbors:'Measure the adjacency statistics for the nuclei. Nuclei within 1 pixel of each other are considered neighbors.'
# # 	Select objects to measure:Nuclei
# #     Select neighboring objects to measure:Nuclei
# #     Method to determine neighbors:Within a specified distance
# #     Neighbor distance:1
# #     Retain the image of objects colored by numbers of neighbors?:No
# #     Name the output image:ObjectNeighborCount
# #     Select colormap:Default
# #     Retain the image of objects colored by percent of touching pixels?:No
# #     Name the output image:PercentTouching
# #     Select a colormap:Default

# # MeasureObjectNeighbors:'Measure the adjacency statistics for the cells. Cells touching each other are considered neighbors.'
# #     Select objects to measure:Cells
# #     Select neighboring objects to measure:Cells
# #     Method to determine neighbors:Adjacent
# #     Neighbor distance:5
# #     Retain the image of objects colored by numbers of neighbors?:No
# #     Name the output image:ObjectNeighborCount
# #     Select colormap:Default
# #     Retain the image of objects colored by percent of touching pixels?:No
# #     Name the output image:PercentTouching
# #     Select a colormap:Default

# appears that CellProfiler calculates FirstClosestDistance, SecondClosestDistance, and AngleBetweenNeighbors
# as closest distance between centers of objects identified as neighbors using distances to perimeter. If no
# neighbor close enough to perimeter, then no distance calculated. Here, I have calculated first_neighbor_distance,
# second_neighbor_distances, and angle_between_neighbors using objects with smallest distance between centers, 
# regardless of distance between perimeters. This produces a single metric for all cells, even if multiple distance 
# thresholds are used to find number of perimeter neighbors.

# these features are dependent on information from the entire field-of-view, thus are not extracted with regionprops

def neighbor_measurements(labeled, distances=[1,10],n_cpu=1):
	from pandas import concat

	dfs = [object_neighbors(labeled,distance=distance).rename(columns=lambda x: x+'_'+str(distance)) for distance in distances]

	dfs.append(closest_objects(labeled,n_cpu=n_cpu).drop(columns=['first_neighbor','second_neighbor']))

	return concat(dfs,axis=1,join='outer').reset_index()

# # MeasureObjectRadialDistribution:'Measure the radial intensity distribution characteristics in all objects. The object is "binned" into radial annuli and statistics are measured for each bin.'
# #     Hidden:4
# #     Hidden:3
# #     Hidden:1
# #     Select an image to measure:ER
# #     Select an image to measure:RNA
# #     Select an image to measure:AGP
# #     Select an image to measure:Mito
# #     Select objects to measure:Cells
# #     Object to use as center?:These objects
# #     Select objects to use as centers:None
# #     Select objects to measure:Nuclei
# #     Object to use as center?:These objects
# #     Select objects to use as centers:None
# #     Select objects to measure:Cytoplasm
# #     Object to use as center?:These objects
# #     Select objects to use as centers:None
# #     Scale the bins?:Yes
# #     Number of bins:4
# #     Maximum radius:100

# This module no longer exists in CellProfiler -> MeasureObjectIntensityDistribution
# But do not use intensity zernike's--slow and not useful: https://github.com/CellProfiler/CellProfiler/issues/2220
# center defined as point farthest from edge = np.argmax(distance_transform(np.pad(r.filled_image,1,'constant')))

intensity_distribution_features = {
	# to minimize re-computing values, outputs a numpy array of length 3*bins. order is [FracAtD, MeanFrac, RadialCV]*bins
	# pretty expensive and dubious utility, leave out if computation is limiting
	'intensity_distribution' : lambda r: np.array(
		measure_intensity_distribution(r.filled_image,r.image,r.intensity_image,bins=4)
		).reshape(-1),
	'weighted_hu_moments': lambda r: catch_runtime(lambda r:r.weighted_moments_hu)(r)
}

intensity_distribution_features_ch = {
	# to minimize re-computing values, outputs a numpy array of length 3*bins. order is [FracAtD, MeanFrac, RadialCV]*bins
	# pretty expensive and dubious utility, leave out if computation is limiting
	'intensity_distribution' : lambda r, ch: np.array(
		measure_intensity_distribution(r.filled_image,r.image,r.intensity_image_full[ch],bins=4)
		).flatten(),
	'weighted_hu_moments': lambda r, ch: weighted_hu_moments_grayscale(r.intensity_image_full[ch]*r.image)
}

intensity_distribution_features_multichannel = {
	# to minimize re-computing values, outputs a numpy array of length 3*bins. order is [FracAtD, MeanFrac, RadialCV]*bins
	# pretty expensive and dubious utility, leave out if computation is limiting
	'intensity_distribution' : lambda r: np.concatenate(
		measure_intensity_distribution_multichannel(r.filled_image,r.image,r.intensity_image,bins=4)),
	'weighted_hu_moments': lambda r: catch_runtime(lambda r: r.weighted_moments_hu)(r).flatten()
}

intensity_distribution_columns = {
	'intensity_distribution_0':'frac_at_d_0',
	'intensity_distribution_1':'frac_at_d_1',
	'intensity_distribution_2':'frac_at_d_2',
	'intensity_distribution_3':'frac_at_d_3',
	'intensity_distribution_4':'mean_frac_0',
	'intensity_distribution_5':'mean_frac_1',
	'intensity_distribution_6':'mean_frac_2',
	'intensity_distribution_7':'mean_frac_3',
	'intensity_distribution_8':'radial_cv_0',
	'intensity_distribution_9':'radial_cv_1',
	'intensity_distribution_10':'radial_cv_2',
	'intensity_distribution_11':'radial_cv_3',
}

intensity_distribution_columns_ch = {
	'{channel}_intensity_distribution_0':'{channel}_frac_at_d_0',
	'{channel}_intensity_distribution_1':'{channel}_frac_at_d_1',
	'{channel}_intensity_distribution_2':'{channel}_frac_at_d_2',
	'{channel}_intensity_distribution_3':'{channel}_frac_at_d_3',
	'{channel}_intensity_distribution_4':'{channel}_mean_frac_0',
	'{channel}_intensity_distribution_5':'{channel}_mean_frac_1',
	'{channel}_intensity_distribution_6':'{channel}_mean_frac_2',
	'{channel}_intensity_distribution_7':'{channel}_mean_frac_3',
	'{channel}_intensity_distribution_8':'{channel}_radial_cv_0',
	'{channel}_intensity_distribution_9':'{channel}_radial_cv_1',
	'{channel}_intensity_distribution_10':'{channel}_radial_cv_2',
	'{channel}_intensity_distribution_11':'{channel}_radial_cv_3',
}

intensity_distribution_columns_multichannel = {
	'intensity_distribution':['frac_at_d_0', 'frac_at_d_1', 
							  'frac_at_d_2', 'frac_at_d_3',
							  'mean_frac_0', 'mean_frac_1',
							  'mean_frac_2', 'mean_frac_3',
							  'radial_cv_0', 'radial_cv_1',
							  'radial_cv_2', 'radial_cv_3'],
	'weighted_hu_moments':[f'weighted_hu_moments_{n}' for n in range(7)]
	}
# # MeasureObjectSizeShape:'Measure the morpholigical features of all objects.'
# #     Select objects to measure:Cells
# #     Select objects to measure:Nuclei
# #     Select objects to measure:Cytoplasm
# #     Calculate the Zernike features?:Yes

shape_features = {
	'area'    : lambda r: r.area,
	'perimeter' : lambda r: r.perimeter,
	'convex_area' : lambda r: r.convex_area,
	'form_factor': lambda r:4*np.pi*r.area/(r.perimeter)**2, #isoperimetric quotient
	'solidity': lambda r: r.solidity,
	'extent': lambda r: r.extent,
	'euler_number': lambda r: r.euler_number,
	'centroid': lambda r: r.local_centroid, # ACTUALLY SHOULD BE POINT FARTHEST FROM EDGE?
	'eccentricity': lambda r: r.eccentricity,
	'major_axis' : lambda r: r.major_axis_length,
    'minor_axis' : lambda r: r.minor_axis_length,
    'orientation' : lambda r: r.orientation,
    'compactness' : lambda r: 2*np.pi*(r.moments_central[0,2]+r.moments_central[2,0])/(r.area**2), # centrosome.cpmorphology.ellipse_from_second_moments(): "variance of the radial distribution normalized by the area"
    'radius' : lambda r: max_median_mean_radius(r.filled_image),
    'feret_diameter' : lambda r: min_max_feret_diameter(r.coords), # relatively expensive, likely high correlation with major/minor axis; looks like max feret will be added to skimage regionprops, but not yet
    'hu_moments': lambda r: r.moments_hu,
    # zernike okay to remove if computation is limiting -- not many of these were retained in Rohban 2017 eLife and they are very (most) expensive
    # cp/centrosome zernike divides zernike magnitudes by minimum enclosing circle magnitude; unclear why
    'zernike' : lambda r: zernike_minimum_enclosing_circle(r.coords, degree=ZERNIKE_DEGREE)
}

zernike_nums = ['zernike_'+str(radial)+'_'+str(azimuthal) 
for radial in range(ZERNIKE_DEGREE+1) 
for azimuthal in range(radial%2,radial+2,2)]

shape_columns = {'zernike_'+str(num):zernike_num for num,zernike_num in enumerate(zernike_nums)}
shape_columns.update({
	'centroid_0':'centroid_r',
	'centroid_1':'centroid_c',
	'radius_0':'max_radius',
	'radius_1':'median_radius',
	'radius_2':'mean_radius',
	'feret_diameter_0':'min_feret_diameter',
	'feret_diameter_1':'max_feret_diameter',
	# these are really for plotting purposes
	'feret_diameter_2':'min_feret_r0',
	'feret_diameter_3':'min_feret_c0',
	'feret_diameter_4':'min_feret_r1',
	'feret_diameter_5':'min_feret_c1',
	'feret_diameter_6':'max_feret_r0',
	'feret_diameter_7':'max_feret_c0',
	'feret_diameter_8':'max_feret_r1',
	'feret_diameter_9':'max_feret_c1'
})

# # MeasureTexture:'Measure the texture features in all objects, against all 5 channels, using multiple spatial scales.'
# #     Hidden:5
# #     Hidden:3
# #     Hidden:3
# #     Select an image to measure:DNA
# #     Select an image to measure:ER
# #     Select an image to measure:RNA
# #     Select an image to measure:AGP
# #     Select an image to measure:Mito
# #     Select objects to measure:Cells
# #     Select objects to measure:Cytoplasm
# #     Select objects to measure:Nuclei
# #     Texture scale to measure:3
# #     Angles to measure:Horizontal
# #     Texture scale to measure:5
# #     Angles to measure:Horizontal
# #     Texture scale to measure:10
# #     Angles to measure:Horizontal
# #     Measure Gabor features?:Yes
# #     Number of angles to compute for Gabor:4

# Gabor features not used in more recent version of cell painting analysis pipeline
# each haralick feature outputs 13 features
# unclear how cell profiler aggregates results from all 4 directions, most likely is mean

# Haralick references:
#	Haralick RM, Shanmugam K, Dinstein I. (1973), “Textural Features for Image Classification” IEEE Transaction on Systems Man, Cybernetics, SMC-3(6):610-621.
#	http://murphylab.web.cmu.edu/publications/boland/boland_node26.html

texture_features = {
	# PFTAS is an alternative to Haralick for texture: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-110
	# mahotas implementation has some quircks, e.g., outputs 54 features for 2D image: the 9 PFTAS statistics for 3 different binary images 
	# and their complement. The first 9 features are the "realest" PFTAS statistics.
	'pftas' : lambda r: masked_pftas(r.intensity_image),
	'haralick_5'  : lambda r: ubyte_haralick(r.intensity_image, ignore_zeros=True, distance=5,  return_mean=True)
	# haralick computational cost increasing significantly with distance; just keep local 5 pixel texture here for most uses.
	# 'haralick_10' : lambda r: ubyte_haralick(r.intensity_image, ignore_zeros=True, distance=10, return_mean=True),
	# 'haralick_20' : lambda r: ubyte_haralick(r.intensity_image, ignore_zeros=True, distance=20, return_mean=True)
}

texture_features_ch = {
	# PFTAS is an alternative to Haralick for texture: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-110
	# mahotas implementation has some quircks, e.g., outputs 54 features for 2D image: the 9 PFTAS statistics for 3 different binary images 
	# and their complement. The first 9 features are the "realest" PFTAS statistics.
	'pftas' : lambda r, ch: masked_pftas(r.intensity_image_full[ch]*r.image),
	'haralick_5'  : lambda r, ch: ubyte_haralick(r.intensity_image_full[ch]*r.image, ignore_zeros=True, distance=5,  return_mean=True)
	# haralick computational cost increasing significantly with distance; just keep local 5 pixel texture here for most uses.
	# 'haralick_10' : lambda r: ubyte_haralick(r.intensity_image, ignore_zeros=True, distance=10, return_mean=True),
	# 'haralick_20' : lambda r: ubyte_haralick(r.intensity_image, ignore_zeros=True, distance=20, return_mean=True)
}

texture_features_multichannel = {
	# PFTAS is an alternative to Haralick for texture: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-110
	# mahotas implementation has some quircks, e.g., outputs 54 features for 2D image: the 9 PFTAS statistics for 3 different binary images 
	# and their complement. The first 9 features are the "realest" PFTAS statistics.
	'pftas' : lambda r: np.array([masked_pftas(channel) 
		for channel in np.moveaxis(r.intensity_image.reshape(*r.intensity_image.shape[:2],-1),-1,0)]).flatten(order='F'),
	# 'haralick_5'  : lambda r: np.array([ubyte_haralick(channel, ignore_zeros=True, distance=5,  return_mean=True)
	# 	for channel in np.moveaxis(r.intensity_image.reshape(*r.intensity_image.shape[:2],-1),-1,0)]).flatten(order='F'),
	# # haralick computational cost increasing significantly with distance; just keep local 5 pixel texture here for most uses.
	'haralick_10'  : lambda r: np.array([ubyte_haralick(channel, ignore_zeros=True, distance=10,  return_mean=True)
		for channel in np.moveaxis(r.intensity_image.reshape(*r.intensity_image.shape[:2],-1),-1,0)]).flatten(order='F'),
	# 'haralick_20'  : lambda r: np.array([ubyte_haralick(channel, ignore_zeros=True, distance=20,  return_mean=True)
	# 	for channel in np.moveaxis(r.intensity_image.reshape(*r.intensity_image.shape[:2],-1),-1,0)]).flatten(order='F')
}

texture_columns_multichannel = {
	'pftas':[f'pftas_{n}' for n in range(54)],
	# 'haralick_5':[f'haralick_5_{n}' for n in range(13)],
	'haralick_10':[f'haralick_10_{n}' for n in range(13)],
	# 'haralick_20':[f'haralick_20_{n}' for n in range(13)]
}

######################################################################################################################################

grayscale_features = {**intensity_features,
					  **intensity_distribution_features,
					  **texture_features,
					  # really slow
					  # **granularity_features
					 }

grayscale_features_ch = {**intensity_features_ch,
					  **intensity_distribution_features_ch,
					  **texture_features_ch,
					 }

grayscale_features_multichannel = {**intensity_features_multichannel,
					  **intensity_distribution_features_multichannel,
					  **texture_features_multichannel,
					 }

grayscale_columns = {**intensity_columns,
					 **intensity_distribution_columns
					}

grayscale_columns_ch = {**intensity_columns_ch,
					 **intensity_distribution_columns_ch
					}

grayscale_columns_multichannel = {**intensity_columns_multichannel,
					**intensity_distribution_columns_multichannel,
					**texture_columns_multichannel
					}
					
######################################################################################################################################

@decorator
def catch_runtime(func,*args,**kwargs):
	with catch_warnings():
		simplefilter("ignore",category=RuntimeWarning)
		return func(*args,**kwargs)

def lstsq_slope(r,first,second):
	A = masked(r,first)
	B = masked(r,second)

	filt = A > 0
	if filt.sum() == 0:
	    return np.nan

	A = A[filt]
	B  = B[filt]
	slope = np.linalg.lstsq(np.vstack([A,np.ones(len(A))]).T,B,rcond=-1)[0][0]

	return slope

def lstsq_slope_all_multichannel(r):
	V = r.intensity_image[r.image]

	slopes = []
	for ch in range(r.intensity_image.shape[-1]):
		slopes.extend(np.linalg.lstsq(np.vstack([V[...,ch],np.ones(V.shape[0])]).T,
			np.delete(V,ch,axis=1),rcond=None)[0][0])

	return slopes

def cp_colocalization_all_channels(r,mode='multichannel',**kwargs):
	if mode=='multichannel':
		channels = r.intensity_image.shape[-1]
	else:
		channels = r.intensity_image_full.shape[-3]

	results = [cp_colocalization(r,first,second,mode=mode,**kwargs) 
				for first,second in combinations(list(range(channels)),2)]

	if mode=='multichannel':
		return np.array(results).flatten(order='F')
	else:
		return np.concatenate(results)

def cp_colocalization(r,first,second,mode='multichannel',**kwargs):
	if mode=='multichannel':
		A, B = r.intensity_image[r.image][...,[first, second]].T
	else:
		A = masked(r,first)
		B = masked(r,second)
	return measure_colocalization(A,B,**kwargs)

def measure_colocalization(A,B,threshold='otsu'):
	"""Measures overlap, k1/k2, manders, and rank weighted colocalization coefficients.
	References:
	http://www.scian.cl/archivos/uploads/1417893511.1674 starting at slide 35
	Singan et al. (2011) "Dual channel rank-based intensity weighting for quantitative 
	co-localization of microscopy images", BMC Bioinformatics, 12:407.
	threshold is either 'otsu' or 'costes' methods, or a float between 0 and 1 defining the 
	fraction of the maximum value to be used as the threshold (CellProfiler default=0.15)
	"""
	if (A.sum()==0) | (B.sum()==0):
		return (np.nan,)*7

	results = []

	if threshold=='otsu':
		A_thresh, B_thresh = otsu(A),otsu(B)
	elif threshold=='costes':
		A_thresh,B_thresh = costes_threshold(A,B)
	elif ifinstance(threshold,float)&(0<=threshold<=1):
		A_thresh,B_thresh = (threshold*A.max(), threshold*B.max())
	else:
		raise ValueError('`threshold` must be a float on the interval [0,1] or one of the methods "otsu" or "costes"')

	A, B = A.astype(float),B.astype(float)

	# commented out versions reflect actual CellProfiler computations, but these don't match literature
	# A_total, B_total = A[A>A_thresh].sum(), B[B>B_thresh].sum()

	# mask = (A > A_thresh) & (B > B_thresh)
	# A_mask = A[mask]
	# B_mask = B[mask] 

	# overlap = (A_mask*B_mask).sum()/np.sqrt((A_mask**2).sum()*(B_mask**2).sum())
	overlap = (A*B).sum()/np.sqrt((A**2).sum()*(B**2).sum())

	results.append(overlap)

	# K1 = (A_mask*B_mask).sum()/(A_mask**2).sum()
	# K2 = (A_mask*B_mask).sum()/(B_mask**2).sum()
	K1 = (A*B).sum()/(A**2).sum()
	K2 = (A*B).sum()/(B**2).sum()

	results.extend([K1,K2])

	# M1 = A_mask.sum()/A_total
	# M2 = B_mask.sum()/B_total
	M1 = A[B>B_thresh].sum()/A.sum()
	M2 = B[A>A_thresh].sum()/B.sum()

	results.extend([M1,M2])

	A_ranks = rankdata(A,method='dense')
	B_ranks = rankdata(B,method='dense')

	R = max([A_ranks.max(),B_ranks.max()])
	# weight = ((R-abs(A_ranks-B_ranks))/R)[mask]
	# RWC1 = (A_mask*weight).sum()/A_total
	# RWC2 = (B_mask*weight).sum()/B_total
	weight = ((R-abs(A_ranks-B_ranks))/R)
	RWC1 = ((A*weight)[B>B_thresh]).sum()/A.sum()
	RWC2 = ((B*weight)[A>A_thresh]).sum()/B.sum()

	results.extend([RWC1,RWC2])

	return results

def costes_threshold(A,B,step=1,pearson_cutoff=0):
	# Costes et al. (2004) Biophysical Journal, 86(6) 3993-4003
	# iteratively decreases threshold until pixels below the threshold have pearson correlation < 0
	# doesn't work if pearson correlation for unthresholded pixels starts is negative
	A_dtype_max, B_dtype_max = np.iinfo(A.dtype).max,np.iinfo(B.dtype).max
	if A_dtype_max != B_dtype_max:
		raise ValueError('inputs must be of the same dtype')
	A = A/A_dtype_max
	B = B/A_dtype_max
	# step = step/xA_dtype_max

	mask = (A>0)|(B>0)

	A = A[mask]
	B = B[mask]

	A_var = np.var(A,ddof=1)
	B_var = np.var(B,ddof=1)

	Z = A+B
	Z_var = np.var(Z,ddof=1)

	covar = 0.5 * (Z_var - (A_var+B_var))

	a = (B_var-A_var)+np.sqrt((B_var-A_var)**2 + 4*(covar**2))/(2*covar)

	b = B.mean()-a*A.mean()

	threshold = A.max()

	if (len(np.unique(A)) > 10**4) & (step<100):
		step = 100

	# could also try the histogram bisection method used in Coloc2
	# https://github.com/fiji/Colocalisation_Analysis
	for threshold in np.unique(A)[::-step]:
		below = (A<threshold)|(B<(a*threshold+b))
		pearson = np.mean((A[below]-A[below].mean())*(B[below]-B[below].mean())/(A[below].std()*B[below].std()))

		if pearson <= pearson_cutoff:
			break

	return threshold*A_dtype_max,(a*threshold+b)*B_dtype_max

# def granularity_spectrum(grayscale, labeled, background_radius=5, spectrum_length=16, downsample=1, background_downsample=0.5):
# 	"""Returns granularity spectrum as defined in the CellProfiler documentation.
# 	Scaled so that units are approximately the % of new granules stuck in imaginary sieve when moving to 
# 	size specified by spectrum component
# 	Helpful resources:
# 	Maragos P. “Pattern spectrum and multiscale shape representation”,
# 		IEEE Transactions on Pattern Analysis and Machine Intelligence, 
# 		VOL 11, NO 7, pp. 701-716, 1989
# 	Vincent L. (1992) “Morphological Area Opening and Closing for
# 		Grayscale Images”, Proc. NATO Shape in Picture Workshop,
# 		Driebergen, The Netherlands, pp. 197-208.
# 	https://en.wikipedia.org/wiki/Granulometry_(morphology)
# 	http://www.ravkin.net/presentations/Statistical%20properties%20of%20algorithms%20for%20analysis%20of%20cell%20images.pdf
# 	"""
# 	intensity_image = grayscale.copy()
# 	image = labeled.copy()


# 	i_sub,j_sub = np.mgrid[0:image.shape[0]*downsample, 0:image.shape[1]*downsample].astype(float)/downsample
# 	if downsample < 1:
# 		intensity_image = map_coordinates(intensity_image,(i_sub,j_sub),order=1)
# 		image = map_coordinates(image.astype(float),(i_sub,j_sub))>0.9

# 	if background_downsample <1:
# 		i_sub_sub,j_sub_sub = (np.mgrid[0:image.shape[0]*background_downsample, 
# 			0:image.shape[1]*background_downsample].astype(float)/background_downsample)
# 		background_intensity = map_coordinates(intensity_image,(i_sub_sub,j_sub_sub),order=1)
# 		background_mask = map_coordinates(image.astype(float),(i_sub_sub,j_sub_sub))>0.9
# 	else:
# 		background_intensity = intensity_image
# 		background_mask = image

# 	selem = skimage.morphology.disk(background_radius,dtype=bool)

# 	# cellprofiler masks before and between erosion/dilation steps here--
# 	# this creates unwanted edge effects here. Combine erosion/dilation into opening
# 	# background = skimage.morphology.erosion(background_intensity*background_mask,selem=selem)
# 	# background = skimage.morphology.dilation(background,selem=selem)
# 	background = skimage.morphology.opening(background_intensity,selem=selem)

# 	# rescaling
# 	if background_downsample < 1:
# 		# rescale background to match intensity_image
# 		i_sub *= float(background.shape[0]-1)/float(image.shape[0]-1)
# 		j_sub *= float(background.shape[1]-1)/float(image.shape[1]-1)
# 		background = map_coordinates(background,(i_sub,j_sub),order=1)

# 	# remove background
# 	intensity_image -= background
# 	intensity_image[intensity_image<0] = 0

# 	# calculate granularity spectrum
# 	start = np.mean(intensity_image[image])

# 	# cellprofiler also does unwanted masking step here
# 	erosion = intensity_image

# 	current = start

# 	footprint = skimage.morphology.disk(1,dtype=bool)

# 	spectrum = []
# 	for _ in range(spectrum_length):
# 		previous = current.copy()
# 		# cellprofiler does unwanted masking step here
# 		erosion = skimage.morphology.erosion(erosion, selem=footprint)
# 		# masking okay here--inhibits bright regions from outside object being propagated into the image
# 		reconstruction = skimage.morphology.reconstruction(erosion*image, intensity_image, selem=footprint)
# 		current = np.mean(reconstruction[image])
# 		spectrum.append((previous - current) * 100 / start)

# 	return spectrum

def boundaries(labeled,connectivity=1,mode='inner',background=0):
    """Supplement skimage.segmentation.find_boundaries to include image edge pixels of 
    labeled regions as boundary
    """
    from skimage.segmentation import find_boundaries
    kwargs = dict(connectivity=connectivity,
        mode=mode,
        background=background
        )
    # if mode == 'inner':
    pad_width = 1
    # else:
    #     pad_width = connectivity

    padded = np.pad(labeled,pad_width=pad_width,mode='constant',constant_values=background)
    return find_boundaries(padded,**kwargs)[...,pad_width:-pad_width,pad_width:-pad_width]

def edge_intensity_features(intensity_image,filled_image,**kwargs):
	edge_pixels = intensity_image[boundaries(filled_image,**kwargs),...]

	return np.array([edge_pixels.sum(axis=0),
		edge_pixels.mean(axis=0),
		np.std(edge_pixels,axis=0),
		edge_pixels.max(axis=0),
		edge_pixels.min(axis=0)
		]
		).flatten()

def weighted_local_centroid_grayscale(intensity_image):
	if intensity_image.sum()==0:
		return (np.nan,)*2
	wm = skimage.measure.moments(intensity_image,order=3)
	return (wm[tuple(np.eye(intensity_image.ndim, dtype=int))] /
                wm[(0,) * intensity_image.ndim])

def weighted_local_centroid_multichannel(r):
	with catch_warnings():
		simplefilter("ignore",category=RuntimeWarning)
		return r.weighted_local_centroid

def mass_displacement_grayscale(local_centroid,intensity_image):
	weighted_local_centroid = weighted_local_centroid_grayscale(intensity_image)
	return np.sqrt(((np.array(local_centroid) - np.array(weighted_local_centroid))**2).sum())

def closest_objects(labeled,n_cpu=1):
	from ops.process import feature_table
	from scipy.spatial import cKDTree

	features = {
	'i'       : lambda r: r.centroid[0],
    'j'       : lambda r: r.centroid[1],
    'label'   : lambda r: r.label
    }
	
	df = feature_table(labeled,labeled,features)

	kdt = cKDTree(df[['i','j']])

	distances,indexes = kdt.query(df[['i','j']],3,workers=n_cpu)

	df['first_neighbor'],df['first_neighbor_distance'] = indexes[:,1],distances[:,1]
	df['second_neighbor'],df['second_neighbor_distance'] = indexes[:,2],distances[:,2]

	first_neighbors = df[['i','j']].values[df['first_neighbor'].values]
	second_neighbors = df[['i','j']].values[df['second_neighbor'].values]

	angles = [angle(v,p0,p1) 
          for v,p0,p1 
          in zip(df[['i','j']].values,first_neighbors,second_neighbors)]

	df['angle_between_neighbors'] = np.array(angles)*(180/np.pi)

	return df.drop(columns=['i','j']).set_index('label')

def object_neighbors(labeled, distance=1):
	from skimage.measure import regionprops
	from pandas import DataFrame
	
	outlined = boundaries(labeled,connectivity=EDGE_CONNECTIVITY,mode='inner')*labeled

	regions = regionprops(labeled)

	bboxes = [r.bbox for r in regions]

	labels = [r.label for r in regions]

	neighbors_disk = skimage.morphology.disk(distance)

	perimeter_disk = cp_disk(distance+0.5)

	info_dicts = [neighbor_info(labeled,outlined,label,bbox,distance,neighbors_disk,perimeter_disk) for label,bbox in zip(labels,bboxes)]

	return DataFrame(info_dicts).set_index('label')

def neighbor_info(labeled,outlined,label,bbox,distance,neighbors_disk=None,perimeter_disk=None):
	if neighbors_disk is None:
		neighbors_disk = skimage.morphology.disk(distance)
	if perimeter_disk is None:
		perimeter_disk = cp_disk(distance+0.5)

	label_mask = subimage(labeled,bbox,pad=distance)
	outline_mask = subimage(outlined,bbox,pad=distance) == label

	dilated = skimage.morphology.binary_dilation(label_mask==label,selem=neighbors_disk)
	neighbors = np.unique(label_mask[dilated])
	neighbors = neighbors[(neighbors!=0)&(neighbors!=label)]
	n_neighbors = len(neighbors)

	dilated_neighbors = skimage.morphology.binary_dilation((label_mask!=label)&(label_mask!=0),selem=perimeter_disk)
	percent_touching = (outline_mask&dilated_neighbors).sum()/outline_mask.sum()

	return {'label':label,'number_neighbors':n_neighbors,'percent_touching':percent_touching}

def cp_disk(radius):
    """Create a disk structuring element for morphological operations
    
    radius - radius of the disk
    """
    iradius = int(radius)
    x, y = np.mgrid[-iradius : iradius + 1, -iradius : iradius + 1]
    radius2 = radius * radius
    strel = np.zeros(x.shape)
    strel[x * x + y * y <= radius2] = 1
    return strel

@catch_runtime
def measure_intensity_distribution(filled_image, image, intensity_image, bins=4):
	if intensity_image.sum()==0:
		return (np.nan,)*12

	binned, center = binned_rings(filled_image,image,bins)

	frac_at_d = np.array([intensity_image[binned==bin_ring].sum() for bin_ring in range(1,bins+1)])/intensity_image[image].sum()

	frac_pixels_at_d = np.array([(binned==bin_ring).sum() for bin_ring in range(1,bins+1)])/image.sum()

	mean_frac = frac_at_d/frac_pixels_at_d

	wedges = radial_wedges(image,center)

	mean_binned_wedges = np.array([np.array([intensity_image[(wedges==wedge)&(binned==bin_ring)].mean() 
		for wedge in range(1,9)]) 
		for bin_ring in range(1,bins+1)])
	radial_cv = np.nanstd(mean_binned_wedges,axis=1)/np.nanmean(mean_binned_wedges,axis=1)

	return frac_at_d,mean_frac,radial_cv

@catch_runtime
def measure_intensity_distribution_multichannel(filled_image, image, intensity_image, bins=4):
	if all((intensity_image[image,...].sum(axis=0))==0):
		return (np.nan,)*12*intensity_image.shape[-1]

	binned, center = binned_rings(filled_image,image,bins)

	frac_at_d = np.array([intensity_image[binned==bin_ring,...].sum(axis=0) for bin_ring in range(1,bins+1)])/intensity_image[image,...].sum(axis=0)

	frac_pixels_at_d = np.array([(binned==bin_ring).sum() for bin_ring in range(1,bins+1)])/image.sum()

	mean_frac = frac_at_d.reshape(bins,-1)/frac_pixels_at_d[:,None]

	wedges = radial_wedges(image,center)

	mean_binned_wedges = np.array([np.array([intensity_image[(wedges==wedge)&(binned==bin_ring),...].mean(axis=0) 
		for wedge in range(1,9)]) 
		for bin_ring in range(1,bins+1)])
	radial_cv = np.nanstd(mean_binned_wedges,axis=1)/np.nanmean(mean_binned_wedges,axis=1)

	return frac_at_d.flatten(),mean_frac.flatten(),radial_cv.flatten()

def binned_rings(filled_image,image,bins):
	"""takes filled image, separates into number of rings specified by bins, 
	with the ring size normalized by the radius at that approximate angle"""

	# returns distance to center point, normalized by distance to edge along
	# that direction, [0,1]; 0 = center point, 1 = points outside the image
	normalized_distance,center = normalized_distance_to_center(filled_image)

	binned = np.ceil(normalized_distance*bins)

	binned[binned==0]=1

	return np.multiply(np.ceil(binned),image),center

def normalized_distance_to_center(filled_image):
	"""regions outside of labeled image have normalized distance of 1"""

	distance_to_edge = distance_transform(np.pad(filled_image,1,'constant'))[1:-1,1:-1]

	max_distance = distance_to_edge.max()

	# median of all points furthest from edge
	center = tuple(np.median(np.where(distance_to_edge==max_distance),axis=1).astype(int))

	mask = np.ones(filled_image.shape)
	mask[center[0],center[1]] = 0

	distance_to_center = distance_transform(mask)

	return distance_to_center/(distance_to_center+distance_to_edge),center

def radial_wedges(image, center):
	"""returns shape divided into 8 radial wedges, each comprising a 45 degree slice
	of the shape from center. Output labeleing convention:
	    i > +
	      \\ 3 || 4 // 
	 +  7  \\  ||  // 8
	 ^  ===============
	 j  5  //  ||  \\ 6
	      // 1 || 2 \\ 
	"""
	i, j = np.mgrid[0 : image.shape[0], 0 : image.shape[1]]

	positive_i,positive_j = (i > center[0], j> center[1])

	abs_i_greater_j = abs(i - center[0]) > abs(j - center[1])

	return ((positive_i + positive_j * 2 + abs_i_greater_j * 4 + 1)*image).astype(int)

def weighted_hu_moments_grayscale(masked_intensity_image):
	if masked_intensity_image.sum()==0:
		return (np.nan,)*7
	return skimage.measure.moments_hu(
		skimage.measure.moments_normalized(
			skimage.measure.moments_central(masked_intensity_image)
			)
		)

def max_median_mean_radius(filled_image):
	transformed = distance_transform(np.pad(filled_image,1,'constant'))[1:-1,1:-1][filled_image]

	return (transformed.max(),
		np.median(transformed),
		transformed.mean()
		)

def min_max_feret_diameter(coords):
	hull_vertices = coords[ConvexHull(coords).vertices]

	antipodes = get_antipodes(hull_vertices)

	point_distances = pdist(hull_vertices)

	try:
		argmin,argmax = (antipodes[:,6].argmin(),point_distances.argmax())
		# min feret diameter, max feret diameter, min feret r0,c0,r1,c1 , max feret r0,c0,r1,c1
		results = ((antipodes[argmin,6],point_distances[argmax])
			+(np.mean([antipodes[argmin,0],antipodes[argmin,2]]),np.mean([antipodes[argmin,1],antipodes[argmin,3]]))
			+tuple(antipodes[argmin,4:6])
			)
		for v in tuple(combinations(hull_vertices,r=2))[argmax]:
			results+=tuple(v)
	except:
		results = (np.nan,)*10

	return results

def get_antipodes(vertices):
    """rotating calipers"""
    antipodes = []
    # iterate through each vertex
    for v_index,vertex in enumerate(vertices):
        current_distance = 0
        candidates = vertices[circular_index(v_index+1,v_index-2,len(vertices))]

        # iterate through each vertex except current and previous
        for c_index,candidate in enumerate(candidates):

            #calculate perpendicular distance from candidate_antipode to line formed by current and previous vertex
            d = perpendicular_distance(vertex,vertices[v_index-1],candidate)
            
            if d < current_distance:
                # previous candidate is a "breaking" antipode
                antipodes.append(np.concatenate([vertex,vertices[v_index-1],candidates[c_index-1],current_distance[None]]))
                break
                
            elif d >= current_distance:
                # not a breaking antipode
                if d == current_distance:
                    # previous candidate is a "non-breaking" antipode
                    antipodes.append(np.concatenate([vertex,vertices[v_index-1],candidates[c_index-1],current_distance[None]]))
                    if c_index == len(candidates)-1:
                        antipodes.append(np.concatenate([vertex,vertices[v_index-1],candidates[c_index],current_distance[None]]))
                current_distance = d

    return np.array(antipodes)

def circular_index(first,last,length):
    if last<first:
        last += length
        return np.arange(first, last+1)%length
    elif last==first:
        return np.roll(range(length),-first)
    else:
        return np.arange(first,last+1)

def perpendicular_distance(line_p0,line_p1,p0):
    if line_p0[0]==line_p1[0]:
        return abs(line_p0[0]-p0[0])
    elif line_p0[1]==line_p1[1]:
        return abs(line_p0[1]-p0[1])
    else:
        return abs(((line_p1[1]-line_p0[1])*(line_p0[0]-p0[0])-(line_p1[0]-line_p0[0])*(line_p0[1]-p0[1]))/
                np.sqrt((line_p1[1]-line_p0[1])**2+(line_p1[0]-line_p0[0])**2))

def zernike_minimum_enclosing_circle(coords,degree=9):
	image, center, diameter = minimum_enclosing_circle_shift(coords)

	return zernike_moments(image, radius=diameter/2, degree=degree, cm=center)

def minimum_enclosing_circle_shift(coords,pad=1):
	diameter,center = minimum_enclosing_circle(coords)

	# diameter = np.ceil(diameter)

	# have to adjust image size to fit minimum enclosing circle
	shift = np.round(diameter/2 - center)
	shifted = np.zeros((int(np.ceil(diameter)+pad),int(np.ceil(diameter)+pad)))
	# shift = np.round(np.array(shifted.shape)/2 - center)
	coords_shifted = (coords + shift).astype(int)
	shifted[coords_shifted[:,0],coords_shifted[:,1]] = 1
	center_shifted = center + shift

	return shifted, center_shifted, np.ceil(diameter)

def minimum_enclosing_circle(coords):
	# http://www.personal.kent.edu/~rmuhamma/Compgeometry/MyCG/CG-Applets/Center/centercli.htm
	# https://www.cs.princeton.edu/courses/archive/spring09/cos226/checklist/circle.html
	hull_vertices = coords[ConvexHull(coords).vertices]

	s0 = hull_vertices[0]
	s1 = hull_vertices[1]

	iterations = 0

	while True:

		remaining = hull_vertices[(hull_vertices!=s0).max(axis=1)&(hull_vertices!=s1).max(axis=1)]

		angles = np.array(list(map(partial(angle,p0=s0,p1=s1),remaining)))
		
		min_angle = angles.min()

		if min_angle >= np.pi/2:
			# circle diameter is s0-s1, center is mean of s0,s1
			diameter = np.sqrt(((s0-s1)**2).sum())
			center = (s0+s1)/2
			break

		vertex = remaining[np.argmin(angles)]

		remaining_angles = np.array(list(starmap(angle,zip([s1,s0],[s0,vertex],[vertex,s1]))))

		if remaining_angles.max() <= np.pi/2:
			# use circumscribing circle of s0,s1,vertex
			diameter,center = circumscribed_circle(s0,s1,vertex)
			break

		keep = [s0,s1][np.argmax(remaining_angles)]

		s0 = keep
		s1 = vertex

		iterations += 1

		if iterations == len(hull_vertices):
			print('maximum_enclosing_circle did not converge')
			diameter = center = None

	return diameter,center

def angle(vertex, p0, p1):
	v0 = p0 - vertex
	v1 = p1 - vertex

	cosine_angle = np.dot(v0, v1) / (np.linalg.norm(v0) * np.linalg.norm(v1))
	return np.arccos(cosine_angle)

def circumscribed_circle(p0,p1,p2):
	# https://en.wikipedia.org/wiki/Circumscribed_circle
	P = np.array([p0,p1,p2])

	Sx = (1/2)*np.linalg.det(np.concatenate([(P**2).sum(axis=1).reshape(3,1),P[:,1].reshape(3,1),np.ones((3,1))],axis=1))
	Sy = (1/2)*np.linalg.det(np.concatenate([P[:,0].reshape(3,1),(P**2).sum(axis=1).reshape(3,1),np.ones((3,1))],axis=1))
	a = np.linalg.det(np.concatenate([P,np.ones((3,1))],axis=1))
	b = np.linalg.det(np.concatenate([P,(P**2).sum(axis=1).reshape(3,1)],axis=1))

	center = np.array([Sx,Sy])/a
	diameter = 2*np.sqrt((b/a) + (np.array([Sx,Sy])**2).sum()/(a**2))
	return diameter,center

def masked_pftas(intensity_image):
	T = otsu(intensity_image,ignore_zeros=True)
	return pftas(intensity_image,T=T)

@catch_runtime
def ubyte_haralick(intensity_image,**kwargs):
	with catch_warnings():
		simplefilter("ignore",category=UserWarning)
		ubyte_image = img_as_ubyte(intensity_image)
	try:
		features = haralick(ubyte_image,**kwargs)
	except ValueError:
		features = [np.nan]*13

	return features
