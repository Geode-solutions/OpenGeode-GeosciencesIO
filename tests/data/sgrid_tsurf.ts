GOCAD SGrid 1 
HEADER {
name: model_A1_box
volume: false
*painted*variable: Density_1e_4
painted: off
*regions*Region_fract*visible: off
*visible*regions: on
regions: on
*regions*Region_fract_center*visible: off
split_faces: false
*volume*solid: true
ascii: off
double_precision_binary: off
}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME Default
PROJECTION Unknown
DATUM Unknown
AXIS_NAME X Y Z
AXIS_UNIT m m m
ZPOSITIVE Elevation
END_ORIGINAL_COORDINATE_SYSTEM
CLASSIFICATION ReservoirGeologicalModel "Geological Models" "Reservoir Model"
AXIS_N 101 101 11 
PROP_ALIGNMENT CELLS
POINTS_OFFSET 0 
POINTS_FILE FRACT_SET_1fractmodel_A1_box__points@@
FLAGS_OFFSET 0 
FLAGS_FILE FRACT_SET_1fractmodel_A1_box__flags@@

REGION Region_fract 0 
REGION Region_fract_center 1 
REGION_FLAGS_ARRAY_LENGTH 112211 
REGION_FLAGS_BIT_LENGTH 3 
REGION_FLAGS_ESIZE 1 
REGION_FLAGS_OFFSET 0 
REGION_FLAGS_FILE FRACT_SET_1fractmodel_A1_box__region_flags@@

PROPERTY 1 Density_1e_4
PROPERTY_CLASS 1 density_1fract
PROPERTY_KIND 1 " Real Number" 
PROPERTY_CLASS_HEADER 1 density_1fract {
kind: Real Number
unit: unitless
pclip: 99
last_selected_folder: Property
low_clip: 0.0001
high_clip: 0.0001
*colormap*standard*colors: 0 0 0 1 85 0 1 1 170 1 1 0 255 1 0 0 
}
PROPERTY_SUBCLASS 1 QUANTITY Float
PROP_ORIGINAL_UNIT 1 none
PROP_UNIT 1 unitless
PROP_NO_DATA_VALUE 1 -99999
PROP_SAMPLE_STATS 1 100000 0.0001 2.43189e-20 0.0001 0.0001
PROP_ESIZE 1 4
PROP_ETYPE 1  IEEE
PROP_ALIGNMENT 1  CELLS
PROP_PAINTED_FLAG_BIT_POS 1  2 
PROP_FORMAT 1 RAW
PROP_OFFSET 1 0
PROP_FILE 1 FRACT_SET_1fractmodel_A1_box_Frac_Density_1e_4@@
END
GOCAD TSurf 1 
HEADER {
*solid*color: 0.196078 0.803922 0.196078 1
visible: on
*isovalues*Version2009*imported_to_version: Version2011_1p1
*isovalues*Version2011*imported_to_version: Version2011_1p1
name: FractureSet_N20_1fract_visualisation
}
GOCAD_ORIGINAL_COORDINATE_SYSTEM
NAME " Pdgm_Epic Local" 
PROJECTION Unknown
DATUM Unknown
AXIS_NAME X Y Z
AXIS_UNIT m m m
ZPOSITIVE Depth
END_ORIGINAL_COORDINATE_SYSTEM
PROPERTY_CLASS_HEADER X {
kind: X
unit: m
}
PROPERTY_CLASS_HEADER Y {
kind: Y
unit: m
}
PROPERTY_CLASS_HEADER Z {
kind: Z
unit: m
is_z: on
}
PROPERTY_CLASS_HEADER vector3d {
kind: Length
unit: m
}
TFACE
VRTX 1 -3940.94287109375 1428.2666015625 1461.34326171875 
VRTX 2 5864.86279296875 -2786.849853515625 -3610.798828125 
VRTX 3 6255.69775390625 -1355.013427734375 1461.34326171875 
VRTX 4 -4331.77783203125 -3.56964111328125 -3610.798828125 
TRGL 1 4 2 
TRGL 1 2 3 
BSTONE 3 
BORDER 5 3 2 
END
