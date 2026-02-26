# Details ----
#' 03c_spots_mods.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 18 Aug 2024
#' 
#' Content: looking individually at elements
#' + Update 24 Feb 2026
#' ++ pulled these models out from the boxplots
#' ++ automated output with a function (silly that I didn't before)
#' -----------

source('helpers/microHelp.R')
library(lme4)
library(emmeans)

# load----
microE <- read_xlsx('data/output/microE.xlsx')
micro <- read_xlsx('data/output/icpms_clean.xlsx')
micro <- micro %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

# and filter----
acE <- microE %>% 
  filter(Species == 'Anoxypristis cuspidata')

pzE <- microE %>% 
  filter(Species == 'Pristis zijsron')

ppE <- microE %>% 
  filter(Species == 'Pristis pristis')

lnrats <- setdiff(ln_rats, "lnNa") # no need for him

# define run----
run_models <- function(data, elements, formula_rhs, 
                       emmeans_vars, mixed = TRUE){
  results <- list()
  
  for(el in elements){
    formula <- as.formula(paste(el, "~", formula_rhs))
    
    model <- if(mixed) lmer(formula, data = data) else glm(formula, data = data)
    
    cat("\n", el, "=====\n")
    print(summary(model))
    
    for(var in emmeans_vars){
      cat("\n", var, "-----\n")
      print(emmeans(model, as.formula(paste("pairwise ~", var))))
    }
    
    results[[el]] <- model
  }
  invisible(results)
}

# lets go!----
## all lmer----

sink("model output/lmer.txt")
run_models(micro, 
           elements = lnrats,
           formula_rhs = "Zone + Region + Species + (1|ID)",
           emmeans_vars = c("Zone", "Region", "Species"))
sink() #now kept as text file!

## glms----
### ac-----
run_models(acE,
           elements = lnrats,
           formula_rhs = "Rsplit", # w maturity?
           emmeans_vars = c("Rsplit"),
           mixed = FALSE)

### pz----
run_models(pzE,
           elements = lnrats,
           formula_rhs = "RivGrp",
           emmeans_vars = c("RivGrp"),
           mixed = FALSE)

### pp----
run_models(ppE,
           elements = lnrats,
           formula_rhs = "RivGrp",
           emmeans_vars = c("RivGrp"),
           mixed = FALSE)


# Manova----
man1 <- manova(cbind(lnLi, lnB, lnMg,	lnMn, lnNi,	lnCu,
                     lnZn,	lnSr, lnBa,	lnPb) ~ Region*Species*Maturity,
               data = microE)
summary(man1) 

man2 <- manova(cbind(lnLi, lnB, lnMg,	lnMn, lnNi,	lnCu,
                     lnZn,	lnSr, lnBa,	lnPb) ~ Region*Maturity + Species,
               data = microE)
summary(man2)


# older code ----

Anova(lm(lnBa ~ Zone + Region + Species*Maturity,
         data = MicroAdj), type = "II")

model <- lm(lnBa ~ Zone + NS + Species*Maturity, data = MicroAdj)
summary(model)

microE %>% 
  group_by(Species, Maturity, Region) %>% 
  summarise(
    meanBa = mean(Ba_Ca, na.rm = TRUE)
  )

jacE <- acE %>%
  filter(Species == "Anoxypristis cuspidata",
         Maturity == 'Immature') %>% 
  mutate(NS = case_when(Region == "NEQ" ~ "N",
                        Region == "CEQ" ~ "S",
                        .default = NS),
         Region = case_when(Region == "WCY" ~ "GoC",
                            Rsplit == "FNEQ" ~ "FNQ",
                            .default = Region))

