models_CV_17 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_17seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_CV.rds")
models_CV_16 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_16seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_CV.rds")
models_CV_15 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_15seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_CV.rds")
models_CV_14 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_14seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_CV.rds")
models_CV_13 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_13seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_CV.rds")

models_fit_15 <- readRDS("C:/Users/cholleem/Documents/ArtificialData_Streambugs_SDMs/B2_outputs_sdms/2000sites_15seed_8pred_0noisetemp_0misdet_NAtoNA_35taxa_4models_/models_FIT.rds")

models_fit_15$GLM$entire_dataset$training$Tipulidae$model$finalModel$converged

models_CV_17$GLM$Split1$training$Tipulidae$model$finalModel$converged
models_CV_17$GLM$Split2$training$Tipulidae$model$finalModel$converged
models_CV_17$GLM$Split3$training$Tipulidae$model$finalModel$converged

models_CV_16$GLM$Split1$training$Tipulidae$model$finalModel$converged

models_CV_15$GLM$Split1$training$Tipulidae$model$finalModel$converged
models_CV_15$GLM$Split2$training$Tipulidae$model$finalModel$converged
models_CV_15$GLM$Split3$training$Tipulidae$model$finalModel$converged

models_CV_15$GAM$Split1$training$Tipulidae$model$finalModel

models_CV_15$GLM$Split1$training$Stratiomyidae$model$finalModel$converged

models_CV_15$GLM$Split1$training$Elmidae$model$finalModel$converged

Stratiomyidae
Elmidae

