# Details ----
#' 03d_spots_vis.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 12 Sep 2024
#' 
#' Content: exploratory pca (plotly)
#' + Update 24 Feb 2026
#' ++ output cleaned up with theme in microHelp.R
#' -----------

source('helpers/microHelp.R')
set_theme(mytheme)
microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

# Cap site----
AC <- microE %>%
  filter(Species == 'Anoxypristis cuspidata') %>% 
  mutate(across(all_of(els), ~ log(. + 1)))

PP <- microE %>%
  filter(Species == 'Pristis pristis') %>% 
  mutate(across(all_of(els), ~ log(. + 1)))

PZ <- microE %>%
  filter(Species == 'Pristis zijsron') %>% 
  mutate(across(all_of(els), ~ log(. + 1)))

## anoxy----
axis = list(showline = FALSE, zeroline = FALSE, 
            gridcolor = '#ffff', ticklen = 5)

fig <- AC %>%  
  plot_ly()  %>%  
  add_trace(  
    type = 'splom',  
    dimensions = list( 
     list(label = 'Li_Ca',values=~Li_Ca),  
     list(label = 'B_Ca',values=~B_Ca),  
     list(label = 'Mn_Ca',values=~Mn_Ca),  
     list(label = 'Mg_Ca',values=~Mg_Ca),  
     list(label = 'Zn_Ca',values=~Zn_Ca),  
     #list(label = 'Cu_Ca',values=~Cu_Ca),  
     #list(label = 'Ni_Ca',values=~Ni_Ca),  
     list(label = 'Sr_Ca',values=~Sr_Ca),  
     list(label = 'Ba_Ca',values=~Ba_Ca)),  
     #list(label = 'Pb_Ca',values=~Pb_Ca),  
     color = ~RivGrp)
     
fig %>%
  layout(
    legend = list(title = list(text = 'Species')),
    hovermode = 'closest', dragmode = 'select',
    plot_bgcolor = 'rgba(240, 240, 240, 0.95)',
    xaxis = axis, yaxis = axis,
    xaxis2 = axis, yaxis2 = axis,
    xaxis3 = axis, yaxis3 = axis,
    xaxis4 = axis, yaxis4 = axis,
    xaxis5 = axis, yaxis5 = axis,
    xaxis6 = axis, yaxis6 = axis,
    xaxis7 = axis, yaxis7 = axis)

## freshie----
axis = list(showline = FALSE, zeroline = FALSE, 
            gridcolor = '#ffff', ticklen = 4)

fig <- PP %>%  
  plot_ly()  %>%  
  add_trace(  
    type = 'splom',  
    dimensions = list( 
      list(label = 'Li_Ca',values=~Li_Ca),  
      list(label = 'B_Ca',values=~B_Ca),  
      list(label = 'Mn_Ca',values=~Mn_Ca),  
      list(label = 'Mg_Ca',values=~Mg_Ca),  
      list(label = 'Zn_Ca',values=~Zn_Ca),  
      #list(label = 'Cu_Ca',values=~Cu_Ca),  
      #list(label = 'Ni_Ca',values=~Ni_Ca),  
      list(label = 'Sr_Ca',values=~Sr_Ca),  
      list(label = 'Ba_Ca',values=~Ba_Ca),  
      list(label = 'Pb_Ca',values=~Pb_Ca)),  
    color = ~ RivGrp)

fig %>%
  layout(
    legend = list(title = list(text = 'Species')),
    hovermode = 'closest', dragmode = 'select',
    plot_bgcolor = 'rgba(240, 240, 240, 0.95)',
    xaxis = axis, yaxis = axis,
    xaxis2 = axis, yaxis2 = axis,
    xaxis3 = axis, yaxis3 = axis,
    xaxis4 = axis, yaxis4 = axis,
    xaxis5 = axis, yaxis5 = axis,
    xaxis6 = axis, yaxis6 = axis,
    xaxis7 = axis, yaxis7 = axis,
    xaxis8 = axis, yaxis8 = axis)

## green----
axis = list(showlinxe=FALSE, zeroline=FALSE, 
            gridcolor='#ffff', ticklen=4)

fig <- PZ %>%  
  plot_ly()  %>%  
  add_trace(  
    type = 'splom',  
    dimensions = list( 
      list(label = 'Li_Ca',values=~Li_Ca),  
      list(label = 'B_Ca',values=~B_Ca),  
      list(label = 'Mn_Ca',values=~Mn_Ca),  
      #list(label = 'Mg_Ca',values=~Mg_Ca),  
      #list(label = 'Zn_Ca',values=~Zn_Ca),  
      #list(label = 'Cu_Ca',values=~Cu_Ca),  
      #list(label = 'Ni_Ca',values=~Ni_Ca),  
      #list(label = 'Sr_Ca',values=~Sr_Ca),  
      list(label = 'Ba_Ca',values=~Ba_Ca),  
      list(label = 'Pb_Ca',values=~Pb_Ca)),  
    color = ~Region)

fig %>%
  layout(
    legend = list(title = list(text = 'Species')),
    hovermode = 'closest', dragmode = 'select',
    plot_bgcolor = 'rgba(240, 240, 240, 0.95)',
    xaxis = axis, yaxis = axis,
    xaxis2 = axis, yaxis2 = axis,
    xaxis3 = axis, yaxis3 = axis,
    xaxis4 = axis, yaxis4 = axis,
    xaxis5 = axis, yaxis5 = axis,
    xaxis6 = axis, yaxis6 = axis,
    xaxis7 = axis, yaxis7 = axis,
    xaxis8 = axis, yaxis8 = axis)
