# Details ----
#' 04a_spots_rf.R
#' Paper: [The potential of elemental analysis in determining origin of sawfish rostra]
#' DOI: tbd
#' Author: Nikki Biskis
#' Date: 5 Sep 2024
#' 
#' Content: random forest models
#' + Update 24 Feb 2026
#' ++ automated for all species
#' +++ was far too long
#' -----------

pacman::p_load('partykit', 'randomForest', 'caret', 
               'rsample', 'rpart', 'e1071')
source('helpers/microHelp.R')

microE <- read_xlsx('data/output/microE.xlsx')
microE <- microE %>% 
  mutate(across(where(is.character) & !all_of(c("ID", "File_ID")), as.factor))

run_rf <- function(data, response, elements, ntree = 500, 
                   mtry = 5, seed = 100, min_n = 3){
  
  dat <- data %>% 
    dplyr::select(all_of(c(response, elements))) %>%
    group_by(.data[[response]]) %>%
    filter(n() >= min_n) %>%
    ungroup() %>%
    droplevels()
  
  set.seed(seed)
  split <- initial_split(dat, prop = 0.7, strata = all_of(response))
  TrainSet <- training(split)
  ValidSet <- testing(split)
  
  model1 <- randomForest(as.formula(paste(response, "~ .")), data = TrainSet, importance = TRUE)
  print(model1)
  
  # or one you set yourself
  model2 <- randomForest(as.formula(paste(response, "~ .")), data = TrainSet, 
                         ntree = ntree, mtry = mtry, importance = TRUE)
  
  predTrain <- predict(model2, TrainSet, type = "class")
  predValid <- predict(model2, ValidSet, type = "class")
  
  cat("\nTrain accuracy:", mean(predTrain == TrainSet[[response]]))
  cat("\nValid accuracy:", mean(predValid == ValidSet[[response]]), "\n")
  print(confusionMatrix(table(predValid, ValidSet[[response]])))
  varImpPlot(model2)
  
  acc <- c()
  for(i in 3:20){
    m <- randomForest(as.formula(paste(response, "~ .")), data = TrainSet, 
                      ntree = ntree, mtry = i, importance = TRUE)
    acc[i-2] <- mean(predict(m, ValidSet, type = "class") == ValidSet[[response]])
  }
  plot(3:20, acc, xlab = "mtry", ylab = "Accuracy", main = paste(response, "- mtry tuning"))
  
  invisible(list(model1 = model1, model2 = model2, valid_acc = acc))
}

# run each species individually, filter/adjust before passing in
AC <- microE %>% filter(Species == "Anoxypristis cuspidata") %>% droplevels()
rf_AC <- run_rf(AC, "Region", els)

PZ <- microE %>% filter(Species == "Pristis zijsron") %>% droplevels()
rf_PZ <- run_rf(PZ, "Region", els, min_n = 2) # too few!

PP <- microE %>% filter(Species == "Pristis pristis") %>% droplevels()
rf_PP <- run_rf(PP, "RivGrp", els)
