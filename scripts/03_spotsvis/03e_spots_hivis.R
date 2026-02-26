# Details ----
#' 03d_spots_vis.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 12 Sep 2024
#' 
#' Content: exploration 2 - plotly tsne and umaps
#' + Update 24 Feb 2026
#' ++ call/output cleaned up with theme in microHelp.R
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))
microH <- microE %>% 
  filter(r <= 100)

# tsne----
library(tsne)

## 2D----
features <- microH[, els]
set.seed(0)
tsne <- tsne(features, initial_dims = 2)
tsne <- data.frame(tsne)

pdb <- cbind(tsne, microH$Species) 
options(warn = -1)
fig <-  plot_ly(data = pdb, x =  ~X1, y = ~X2, 
                type = 'scatter', mode = 'markers', 
                split = ~microH$Species)
fig %>%
  layout(
    plot_bgcolor = "white"
  )

### add Region----
pdb <- cbind(tsne, microH$Species, microH$Region)

colnames(pdb) <- c("X1", "X2", "Species", "Region")
plot_ly(data = pdb, x = ~X1, y = ~X2, 
               color = ~Species, colors = spec_col,
               symbol = ~Region, symbols = 1:13) %>%
  add_markers() %>%
  layout(
    plot_bgcolor = "white",
    legend = list(title = list(text = 'Species')),
    xaxis = list(title = "t-SNE Dimension 1"),
    yaxis = list(title = "t-SNE Dimension 2")
  )

## 3D----
tsne <- tsne(features, initial_dims = 3, k =3)
tsne <- data.frame(tsne)

pdb <- cbind(tsne, microH$Species, microH$ID)
options(warn = -1)
plot_ly(data = pdb ,x =  ~X1, y = ~X2, z = ~X3, 
                color = ~microH$Species, colors = spec_col, text = ~microH$ID) %>% 
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene = list(bgcolor = "white"))

### Region----
pdb <- cbind(tsne, microH$Species, microH$Rsplit, microH$ID)
colnames(pdb) <- c("X1", "X2", "X3", "Species", "Rsplit", "ID")

options(warn = -1)
plot_ly(data = pdb, x =  ~X1, y = ~X2, z = ~X3, 
                symbol = ~Species, color = ~Rsplit, text = ~ID,
        symbols = c('circle', 'square', 'diamond')) %>% 
  add_markers(size = 8) %>%
  layout( 
    xaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'), 
    yaxis = list(
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    scene =list(bgcolor = "white"))

# UMAP----
library(umap) 

umapdat <- microE[, grep("Li_Ca|B_Ca|Mn_Ca|Mg_Ca|Zn_Ca|Cu_Ca|Ni_Ca|Sr_Ca|Ba_Ca|Pb_Ca", 
                          colnames(microE))]
labels <- microE[, "Species"]

spec.umap = umap(umapdat, n_components = 2, random_state = 15) 
layout <- spec.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, microE$Species) 

plot_ly(final, x = ~X1, y = ~X2, color = ~microE$Species, 
               colors = spec_col, type = 'scatter', mode = 'markers')%>%  
  layout(plot_bgcolor = "white",
    legend= list(title=list(text='Species')), 
    xaxis = list(title = "0"), yaxis = list(title = "1")
    ) 

## with region----
spec.umap = umap(umapdat, n_components = 3, random_state = 15) 
layout <- spec.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, microE$Species, microE$Rsplit) 

fig2 <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, 
                symbol = ~microE$Species, color = ~microE$Rsplit,
                symbols = c('circle', 'square', 'diamond'))
fig2 <- fig2 %>% add_markers() 

fig2 %>% layout(scene = list(xaxis = list(title = '0'),
                             yaxis = list(title = '1'), 
                             zaxis = list(title = '2'))) 

 
