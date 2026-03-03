GENERAL INFORMATION
This README.txt file was updated on  9/24/24

A. Paper associated with this archive
Citation: Craig R. McClain, Granger W. Hanks, Mark Genung, Philip J. Manlick, S. River D. Bryant, Avery S. Hiley, James R. Junker, Chiara Romano, Greg W. Rouse, John P. Whiteman, Seth D. Newsome (2024) Abundance, Niche Diversity, and Differential Species Effects Impact Ecosystem Function on the Bottom of the Ocean. American Naturalist
Brief abstract: Biodiversity-ecosystem function (BEF) research stems back decades and has seen a recent exponential growth of the field. Many questions remain about the importance of functional diversity (complementarity) versus redundancy, differential species effects, and the effects of abundance versus different diversity components.  We examined the BEF relationship in deep-sea xylophagous bivalve communities using 63 experimental wood falls deployed ~2000m deep in the Gulf of Mexico. We quantified the relationships between spatial and trophic partitioning, species dominance, individual species effects, and community level properties on total wood consumption, our metric for ecosystem function. A total of 26,324 individuals from 12 xylophagous bivalve species were identified. Individual species inhabited complementary spatial and trophic niches, but species effects on total wood consumption greatly varied. The amount of wood consumed increased significantly with total abundance and species richness, although the richness effects reflected increases in abundance. Ecosystem function among wood fall appeared to be predominately a feature of highly abundant core species with greater functional capacity.

B. Originators
To be included after acceptance

C. Contact information
To be included after acceptance

D. Dates of data collection
In May 2017, 63 wood falls were deployed at two sites (Site 1: N = 23, 2171m and Site 2: N = 40m 1970m) in the northern Gulf of Mexico. The first set of wood falls (Site 1: N=23 and Site 2: N=20) were collected in February/April of 2019 after 21/23 months. The remaining 20 wood falls at Site 2 were collected in February of 2020 after 33 months. 

E. Geographic Location(s) of data collection
Site WF1
27.135N
89.924W
2171m
 
Site WF4
28.814N
87.820W
1970m

F. Funding Sources
NSF: Biological Oceanography.

ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
CC BY-NC 4.0 DEED
2. Data derived from other sources
None
3. Recommended citation for this data/code archive
To be included after acceptance

DATA & CODE FILE OVERVIEW
This data repository consist of 4 data files, 1 code script, and this README document, with the following data and code filenames and variables

Data files and variables

File: woodfall_logs_ef.csv
Description: Data sheet provides information for individual wood falls including size and deployment information
rows: individual wood falls 
column descriptions:
tag number: Number identifier on the plastic tag on each wood fall
tag color: Color of plastic tag on each wood fall
experiment: Type of experiment the wood fall was deployed for.  Site 1-Site 4 refer to the experiments testing 12 woodfalls each of Pinus elliotti and Celtis laevigata varying from small to large total weights. aor refers to experiments testing abundance occupancy relationships (two deployments of 12 small and 12 large wood falls).  wood_type refers to an experiment testing 11 species of wood species on community composition with each species represented by four wood fall sizes.  net_nonet referes to an experiment the effects of the mesh wrap and lack of mesh wrap on the wood falls.  surf_area deployes woodfalls with wood-fall weight controlled but of varying surface area. shallow is used to designate any wood fall deployed 50m and shallower.  
rov dive: Number and letter designation of ROV and dive number for deployment.  GE refers to Oceaneering's Global Explorer.  NA designates no ROV used.
site: Given name for each deployment site.  WF stands for wood fall.
latitude: Deployment site in decimal degrees N˚ 
longitude: Deployment site in decimal degrees W˚
depth: Deployment site in meters
deployment row and column: Arrangement of wood fall deployments at site.  Individual wood falls are deployed 3-4 meters apart. 
wood type: Wood species of each wood fall.  
length: Length of wood fall in centimeters
max diameter: Of the two cut ends, the diameter measurement of the largest in centimeters
min diameter: Of the two cut ends, the diameter measurement of the smallest in centimeters
wood mass final: Final deployed weight of the wood fall. Does not include mesh, zip ties, tags, or rope handles. In kilograms.
notes: Denotes either wood fall categorical size for wood type experiment or aor.  Indicates absence or presence of mesh for net_nonet.  Indicates 2x12 inch plank or 6x6 inch post for surf_area.
weight_wood_zip: Weight of wood fall with zip tie bands. In kilograms. Precision to 0.01 kg. 
weight_packet_total_1: Weight of wood fall with zip tie bands plus mesh covering, additional zip ties, rope handle, and tags. In kilograms. Precision to 0.01 kg. 
weight_wood_only:  Weight of wood fall only with subtracted bands, mesh, zip ties, rope handle, and tags. In kilograms. Precision to 0.01 kg. 
weight_wrap_only: Weight of bands, mesh, zip ties, rope handle, and tags without wood. In kilograms. Precision to 0.01 kg. 
weight_zip: Total weight of zip ties only calculated from standard set of weighed zip ties of zip tie length. In kilograms. Precision to 0.001 kg. 
zip tie length: Length category of zip tie bands used on wood falls. Two zip tie bands per wood fall.  A list of two numbers for a wood fall designates to different sized zip tie bands.
Zip Tie Code	Length (inches)	Mass (kg)
1	60	0.0.255
2	48	0.0214
3	41	0.0192
4	18	0.0077
5	15.5	0.0059
9	29.5	0.01
weight_packet_total_2:  Weight of wood fall with zip tie bands plus mesh covering, additional zip ties, rope handle, and tags. In kilograms.  Weight taken on 5/14/17. Precision to 0.01 kg.

File: woodfall_taxa_ef.csv
rows: individual species
column descriptions:
Description: Data sheet provides taxonomic information for species on wood falls
Phylum/Class/Order/Family/Genus/Species: Updated taxonomic information for each species.  Based on World Registry of Marine Species taxonomy.

File: woodfall_logbytaxa_ef.csv
Description: Data sheet provides counts of species per logs.  Wood fall codes (rows) and taxonomy codes (columns) coincide with codes in woodfall_taxa_ef.csv and woodfall_logs_ef.csv

File: woodfall_taxa_masses_ef.csv
Description: Contains weights of (biomass in mg) of xylophagous bivalves.  Species ID coincides (rows) with species listed in woodfall_taxa_ef.csv.  
Columns: Wood fall number, tag color, and site coincide with woodfall_logs_ef.csv.  Number of individual weighed in biomass quantification is included.  Notes includes which individuals where also used for stable isotope and genetic analysis.

File: WF1_L-G-19_Quarter1_Coordinates.csv
Description: Data sheet provides spatial coordinates of xylophagid bivalves within a single wood fall.  
Rows: single individuals of xylophagid bivalves within the wood fall
Columns: Xylo_group is taxonomy codes (columns) and coincide with codes in woodfall_taxa_ef.csv and woodfall_logs_ef.csv. X is medial distance with larger values indicating deeper burrows, Y is proximal distance from the end of the wood fall along its long axis with smaller values indicating preference for wood fall ends, Z is measured as distance from a tangent line taken from the side of the wood fall used as a metric of surface preference with higher values indicating preference for the top of the wood fall, and W is the maximum burrow width. Y2 is Y recoded as distance to the closer end of the wood fall.

File: XyloTeredi.nex 
Nexus file for xylophagid bivalve phylogeny.

File: wdf_xylos.csv
Description: Data sheet provides isotopic data for xylophagid bivalves in study.
Rows: Individuals on of xylophagid bivalves 
Columns: 
Unique_ID: Combination ID including wood fall number, species ID, and individual count from that wood fall.  Wood fall number and species ID coincide with codes in woodfall_taxa_ef.csv and woodfall_logs_ef.csv
ID: Combination ID including wood fall number, species ID, and individual count form that wood fall.
N.ind: Number of individuals sampled.
species: ID of species, coincide with codes in woodfall_taxa_ef.csv 
log: ID of wood fall, coincide with codes in woodfall_logs_ef.csv
type: Species of tree for wood fall
density: Whether the species of tree is soft or hard wood.
d15N, d13C, percentN. percentC, C.N: isotopic values for individuals.
richness: richness of species on the wood fall
sampled: number of species sampled from the wood fall.

File: wdf_wood.csv
Description: Data sheet provides isotopic and compositional data for wood samples from deployed wood falls in the study.
 Rows: Individual subsamples from each wood fall (from core, mid, and outer sections of each wood type).
 Columns:
	•	Sample_ID: Unique identifier for each subsample, including wood fall ID, tree species code, and core/mid/outer location.
	•	Location: Section of the wood sampled (core, mid, outer).
	•	Species: Scientific name of the tree species.
	•	name: Common name of the tree species.
	•	type: Wood type (hard or soft).
	•	subsample: Type of subsample (sawdust).
	•	d13C: Carbon isotopic value (‰).
	•	%C: Percent carbon.
	•	d15N: Nitrogen isotopic value (‰).
	•	%N: Percent nitrogen.
	•	C_N: Carbon-to-nitrogen ratio.

Code scripts and workflow

File: final.ecosystemfunction2.R 
Description:Complete code for non-isotopic analyses and plots in the paper.

File: WDF_xylos_NicheMetrics.R 
Description: Complete code for isotopic analyses and plots in the paper

File: Niche_functions.R
Description: Functions needed by WDF_xylos_NicheMetrics.R to conduce analyses

File: Niche_functions.R
Description: Complete code for isotopic null model analyses and plots in paper and supplemental

SOFTWARE VERSIONS
R version 4.2.3 (2023-03-15) 
Necessary packages:
Packages (version & purpose):
	•	dplyr 1.0.10 – Data wrangling and manipulation
	•	tidyr 1.2.1 – Reshaping and tidying data
	•	stringr 1.5.0 – String manipulation and pattern matching
	•	tibble 3.1.8 – Enhanced data frames
	•	data.table 1.14.6 – High-performance data manipulation
	•	plyr 1.8.8 – Utility functions for data aggregation
	•	Hmisc 4.8‑1 – Summary statistics and utility functions
	•	ggplot2 3.4.2 – Core plotting and visualization
	•	ggrepel 0.9.2 – Avoid overlapping labels in plots
	•	gridExtra 2.3 – Arrange multiple plots in a grid
	•	viridis 0.6.2 – Color palettes for plots
	•	ggpubr 0.6.0 – Publication-ready ggplot functions
	•	ggthemes 4.2.4 – Additional ggplot themes
	•	vegan 2.6‑4 – Community ecology, diversity metrics
	•	MuMIn 1.46.0 – Model selection and multi-model inference
	•	nlme 3.1‑162 – Linear and nonlinear mixed-effects models
	•	car 3.1‑2 – Companion functions for regression diagnostics
	•	performance 0.9.1 – Model checking and diagnostics
	•	emmeans 1.8.6 – Estimated marginal means (post-hoc tests)
	•	ggeffects 1.1.3 – Visualizing model predictions
	•	candisc 0.4‑8 – Canonical discriminant analysis
	•	mgcv 1.8‑41 – Generalized additive models (GAMs)
	•	SIBER 2.1 – Stable isotope Bayesian ellipses and niche metrics (SEAc, TNW)
	•	spatstat 3.3‑0 – Spatial analysis, polygon/ellipse area calculations


