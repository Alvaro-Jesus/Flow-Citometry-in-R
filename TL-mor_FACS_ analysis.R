#----------------------------------------
#-1. 
install.packages("BiocManager")
# Install cytoinstaller
remotes::install_github("RGLab/cytoinstaller")
# Install cytoverse packages
cytoinstaller::install_cyto(bioc_ver = "devel")
#----------------------------------------


# CytoExploreRData 
devtools::install_github("DillonHammill/CytoExploreRData")
# CytoExploreR 
devtools::install_github("DillonHammill/CytoExploreR")


################### Libraries ###################
pacman::p_load(
  cytolib, 
  flowCore, 
  flowWorkspace, 
  openCyto, 
  ggridges, 
  ggh4x, 
  CytoExploreR,
  tidyverse, 
  tidyplots,
  doBy,
  flowAI,
  shiny,
  ggridges,
  RColorBrewer,
  ggcyto,
  flowAI,
  gridExtra,
  plotly
)



################### Directoy #################


dir_path = "C:/Users/ATL2024_03/OneDrive - Kagoshima University (1)/Maestría/Tesis 2.0/Resultados_Experimentos/Artificial Activation/FACS/TL-mor"
getwd()
setwd(dir_path)

cf_path = "C:/Users/ATL2024_03/OneDrive - Kagoshima University (1)/Maestría/Tesis 2.0/Resultados_Experimentos/Artificial Activation/FACS/TL-mor/Raw"


################## Checking Gating and Log Transformation #########################

#Primitive GatingSet
fs1<-cyto_setup(cf_path, clean=T, sample=T)
fs1
cyto_plot(fs1[1], parent="root", channels = c("FSC-A", "SSC-A"), xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_plot(fs1[4], parent="root", channels = c("FSC-A", "SSC-A"), xlim=c(0, 1e6), ylim=c(0, 1e6))

sampleNames(fs1)

### 4.1 Apply transformation

#Transformed GatingSet
fs1_trans <- cyto_transformer_logicle(fs1, m=5)
fs1t <- cyto_transform(fs1, fs1_trans, plot=F)
sampleNames(fs1t)
class(fs1t)


################ Gating Properly ############################


########### CD25 ###########################

CD25 = cyto_select(fs1t, c(1,2))
sampleNames(CD25)
cyto_gate_draw(CD25, parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD25, parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD25, parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")


#NC_cd25.fcs
cyto_gate_draw(CD25[1], parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD25[1], parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD25[1], parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD25[1], parent="viable", alias="CD25", channels = c("FSC-A", "B585-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")

#PC_cd25.fcs
cyto_gate_draw(CD25[2], parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD25[2], parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD25[2], parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD25[2], parent="viable", alias="CD25", channels = c("FSC-A", "B585-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")




cyto_plot_save("Plot.png",
               height = 7,
               width = 7,
               units = "in",
               res = 300)
cyto_plot(
  x = CD25,
  parent = "viable",
  channels = "B585-A",
  plot_type = "density",
  xlim=c(10^4, 10^7),
  axes_limits = "auto",
  density_stack = 0.15,
  xlab= "CD25 B585-A",
  ylab= NA,
  density_fill_alpha = .4,
  density_cols = c( "#2c7fb8", "gray"),
  title= "TL-mor (HTLV-1 infected cell line)",
  label_text = c("PMA untreated", "PMA treated"),
  density_line_col = c( "black", "black"),
  label_fill_alpha = 0,
  axes_text = c(TRUE, FALSE),
  title_text_size = 2,
  label_text_size = 2,
  axes_text_size = 1.8,
  
  
)



#720 590
############# CD69 #####################


CD69 = cyto_select(fs1t, c(3,4))
sampleNames(CD69)
gs_get_pop_paths(CD69)

cyto_gate_draw(CD69, parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD69, parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD69, parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")


cyto_gate_draw(CD69[1], parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD69[1], parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD69[1], parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD69[1], parent="viable", alias="CD69", channels = c("FSC-A", "B585-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")


cyto_gate_draw(CD69[2], parent="root", alias="cells", channels = c("FSC-A", "SSC-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv", xlim=c(0, 1e6), ylim=c(0, 1e6))
cyto_gate_draw(CD69[2], parent="cells", alias="singlets", channels = c("FSC-A", "FSC-H"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD69[2], parent="singlets", alias="viable", channels = c("FSC-A", "R660-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")
cyto_gate_draw(CD69[2], parent="viable", alias="CD69", channels = c("FSC-A", "B585-A"), axes_limits = "auto", gatingTemplate = "gatingtemp.csv")


#Population moviong = Calculate MFI


cyto_plot_save("Plot2.png",
               height = 7,
               width = 7,
               units = "in",
               res = 300)


cyto_plot(
  x = CD69,
  parent = "viable",
  channels = "B585-A",
  plot_type = "density",
  xlim=c(10^3, 10^6),
  axes_limits = "auto",
  density_stack = 0.15,
  xlab= "CD69 B585-A",
  ylab= NA,
  density_fill_alpha = .4,
  density_cols = c( "#2c7fb8", "gray"),
  title= "TL-mor (HTLV-1 infected cell line)",
  label_text = c("PMA untreated", "    PMA treated"),
  density_line_col = c( "black", "black"),
  label_fill_alpha = 0,
  axes_text = c(TRUE, FALSE),
  title_text_size = 2,
  label_text_size = 2,
  axes_text_size = 1.8,
  
  
)

#720 590

