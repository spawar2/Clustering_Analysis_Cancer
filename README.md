# Clustering_Analysis_Cancer.
[R: GEOquery, Biobase, preprocessCore, multiClust, ctc, gplots, dendextend, graphics, grDevices, amap].
[Paper Link:](https://link.springer.com/chapter/10.1007%2F978-3-030-45385-5_52)
[Poster Link:](https://www.claflin-computation.com/lab-journey?pgid=ktmii98q-ad0f6a9d-26d0-4daf-8ec1-ca96f20fba3c)
<img width="209" alt="Screenshot 2023-05-17 at 10 27 42 AM" src="https://github.com/spawar2/Clustering_Analysis_Cancer/assets/25118302/c212cdd7-4fb1-4f41-80a9-2c9d7cd31a01">
<img width="844" alt="IWBBIO" src="https://github.com/spawar2/Clustering_Analysis_Cancer/assets/25118302/a8420e56-43b9-463a-943b-b6585662a13d">
<img width="159" alt="Clustering" src="https://github.com/spawar2/Clustering_Analysis_Cancer/assets/25118302/8ce6753f-7fe8-4e3d-b4d4-4111a071649b">

[8th International Work-Conference on Bioinformatics and Biomedical Engineering, Granada, Spain, 2020. Clustering reveals common check-point and growth factor receptor genes expressed in six different cancer types, by Pawar SD, Stanam A and Lahiri C.](https://iwbbio.ugr.es/).
Springer Bioinformatics and Biomedical Engineering.
https://campuspress.yale.edu/shrikantpawar/files/2024/05/8-IWBBIO.pptx
https://www.youtube.com/watch?v=Y6skvhHVR2w&ab_channel=ShrikantPawar
Claflin University, Orangeburg, South Carolina, USA. 
https://www.claflin.edu/
https://www.claflin.edu/academics-research/schools-departments/school-of-natural-sciences-and-mathematics/department-of-mathematics-computer-science/computer-science

Cancer_Clustering.R: Breast, Colon, Lung, Oesophageal, Multiple Myeloma, Ovarian Microarray data read, robust multi array (RMA) Normalization, Kmeans analysis, Hierarchal clustering, Plotting.
function(merge, cluster_analysis, hclust, cutree, rbind, heatmap.2, setwd, read.csv, library, set.seed, sample.split, subset, na.omit, scale, svm, predict, table, plot)
