The codes and models for DeepLearning-based algorithm.
# The description of each source code

### Subtyping
This folder contains data and files for subtyping in manuscript. The unet_sample.py and classify.py are used for cell nuclei and area recognition and classification. And xxx.py were further utilized for immune infiltration analysis. The DL_Subtyping.py was used for generating consistant molecular subtyping from proteomic, phosphoproteomic, transciptomic and clinical data.

### Marker
This folder contains data and files for biomarker prediction in manuscripts. The DL_Marker.py was used for generating biomarker combinations for SI, SII and SIII subtypes.

### Kinase
This folder contains data and files for actionable kinase prediction in manuscripts. The DL_Kinase.py was used for generating tumor-specific kinases and subtype-specific kinases.

### SCSP
This folder contains data and files for single cell spatial proteome analysis in manuscripts. The seurat.R was used for generating the initial clusters from cell nuclei image features. The NMF_final.py was used for generating the profile of SCSP. 

# Computation Requirements
### OS Requirements
Above codes have been tested on the following systems: <br>
Windows: Windos 10<br>
Linux: CentOS linux 7.8.2003

### Software Requirements
R (version 4.0.3 or later) tool with packages (Seurat-4.3.0), Python (version 3.8 or later) with modules (Keras-2.4.3, Tensorflow-2.3.1, sklearn-1.1.1).

### Hardware Requirements
All codes and softwares could run on a "normal" desktop computer, no non-standard hardware is needed.<br>
<br>
