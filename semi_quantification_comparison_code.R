source("functions/compound_eluent.R")

# #descriptors and the regressor from previously developed model
#--------------reading in the model-----
descs_pos <-  read_rds("ESIpos_model_descs_191116.rds")
regressor_pos <- read_rds("ESIpos_model_191116.rds")

#----------------graph parameters---------------------------
font <- "Arial"
fontsize <- 9
fontcolour <- "#21313B"
linecolour <- "#21313B"

theme_light <-   theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size = 0.5, colour = fontcolour),
  plot.title = element_text(color = fontcolour,
                            size = 14,
                            face = "bold"),
  text = element_text(family = font,
                      size = fontsize,
                      colour = fontcolour),
  # legend.text = element_text(family = font,
  #                            size = fontsize,
  #                            colour = fontcolour),
  legend.key = element_blank(),
  legend.position = "none",
  axis.text = element_text(family = font,
                           size = fontsize,
                           colour = fontcolour),
  legend.title = element_blank(),
  strip.text = element_text(family = font,
                            size = fontsize),
  strip.background = element_rect(colour = fontcolour,
                                  fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  aspect.ratio = 1,
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))) 

theme_light_wide <-   theme(
  plot.background = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(size = 0.5, colour = fontcolour),
  plot.title = element_text(color = fontcolour,
                            size = 14,
                            face = "bold"),
  text = element_text(family = font,
                      size = fontsize,
                      colour = fontcolour),
  # legend.text = element_text(family = font,
  #                            size = fontsize,
  #                            colour = fontcolour),
  legend.key = element_blank(),
  legend.position = "none",
  axis.text = element_text(family = font,
                           size = fontsize,
                           colour = fontcolour),
  legend.title = element_blank(),
  strip.text = element_text(family = font,
                            size = fontsize),
  strip.background = element_rect(colour = fontcolour,
                                  fill = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks = element_blank(),
  aspect.ratio = 0.8,
  axis.title.x = element_text(hjust = c(1), vjust = c(0)),
  axis.title.y = element_text(hjust = c(1), vjust = c(1))) 

#----------------reading in data------------------------------------
#dataset
dataset <-  read_delim('Dataset_micropollutants.csv',
                          delim = ",",
                          col_names = TRUE,
                          trim_ws = TRUE)

#read the desriptors
Padel_data <-  read_delim('descs.csv',
                          delim = ",",
                          col_names = TRUE,
                          trim_ws = TRUE) %>%
  na.omit()

dataset <- dataset %>%
  left_join(Padel_data) 

#reading the gradient program for reversed phase LC
eluent_parameters <- read_delim('eluent.csv',
                                delim = ",",
                                col_names = TRUE)

#organic modifier in LC gradient
organic_modifier <- "MeOH"
#pH of the waterphase
pH <- 2.7
#NH4 content yes/no
NH4 <- 0

parent <- read_delim('Parents.csv',
                     delim = ",",
                     col_names = TRUE) %>%
  na.omit()

standards <- read_delim('SelectionCalibrationCompounds.csv',
                        delim = ',',
                        col_names = TRUE)

dataset <- dataset %>%
  filter(Sample_type == "Cal Std") %>%
  group_by(Compound) %>%
  mutate(ISTD_area_mean = mean(ISTD_response)) %>%
  ungroup() %>%
  mutate(Signal_corr = Signal/(ISTD_response/ISTD_area_mean)) %>% #correct signal with the ISTD signal
  select(Compound, SMILES, Signal_corr, everything())

#----------combine all data together-------

dataset <- dataset %>%
  mutate(
    c_M = Concentration* (10^-9) /MW, #the starting concentration is in ng/l
    Signal_corr = Signal_corr*IC, #here we account for the isotope peaks
    RF = Signal_corr/c_M,
    organic_modifier = organic_modifier,
    organic = organicpercentage(eluent_parameters, RT),
    pH.aq. = pH,
    NH4 = NH4,
    viscosity =  viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) %>%
  #though all columns are selected with "everything()" in the end, the order of the columns is changed
  dplyr::select(Compound, Type, Sample_type, SMILES, RF,c_M, Signal_corr, RT, MW, IC, organic,viscosity,surface_tension,polarity_index,everything())

#------predicting IEs-------------
prediction_set_model_pos <- dataset %>%
  na.omit() %>%
  mutate(logIE_pred = 0)
prediction <-  predict(regressor_pos, newdata = prediction_set_model_pos, predict.all = TRUE)
prediction <- prediction$aggregate
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  select(SMILES,logIE_pred, everything())

#average the predicted IE values and measured RF values, for calibration solutions.
dataset_summary <- prediction_set_model_pos %>%
  filter(Sample_type == "Cal Std") %>%
  group_by(Compound, SMILES, Type, ISTD_area_mean) %>%
  summarise(logIE_pred = mean(logIE_pred),
            RF = mean(RF),
            RT = mean(RT)) %>%
  ungroup()

#--------choosing the calibration compounds ---------
standards <- standards %>%
  left_join(dataset_summary)

set.seed(0)
split <- sample.split(standards$Compound,
                      SplitRatio = 20) #choose 20 compounds for calibration

#add the calibration compounds information to the matrix
standards <- standards %>%
  mutate(SPLIT = split)

dataset_summary <- dataset_summary %>%
  left_join(standards %>% select(Compound, SPLIT)) %>%
  replace_na(list(SPLIT = FALSE))

cal_compounds <- dataset_summary %>% 
  filter(SPLIT == TRUE)

#linear regression for the calibration compounds
lin_fit_logRF <- lm(log(RF, 10) ~ logIE_pred, 
                    data = dataset_summary %>% filter(SPLIT == TRUE))

dataset_summary <- dataset_summary %>%
  mutate(logRF_pred = lin_fit_logRF$coefficients[2]*logIE_pred + lin_fit_logRF$coefficients[1])
#vis----
#RF vs predicted ionization efficiency
#the calibration compounds are with yellow
ggplot() +
  geom_point(data = dataset_summary %>% filter(SPLIT == FALSE),
             mapping = aes(x = 10^logIE_pred, y = RF), 
             color = "#21313B",
             size = 3, 
             alpha = 1/5) +
  geom_point(data = dataset_summary %>% filter(SPLIT == TRUE),
             mapping = aes(x = 10^logIE_pred, y = RF), 
             color = "#FBC438",
             size = 3, 
             alpha = 4/5) +
  scale_x_log10(limit = c(1e1, 1e5)) +
  scale_y_log10(limit = c(5e17, 5e21)) +
  geom_abline(slope = lin_fit_logRF$coefficients[2], intercept = lin_fit_logRF$coefficients[1], color = "#FBC438", size = 1) +
  scale_shape(solid = FALSE) +
  labs(y = "Response factor", x = "Predicted IE") +
  theme_light +
  annotation_logticks(colour = "gray")

#RF vs RT of the calibration compounds
ggplot() +
  geom_line(data = dataset_summary %>% filter(SPLIT == TRUE),
             mapping = aes(x = RT, y = RF), 
             color = "#21313B",
             size = 1, 
             alpha = 4/5) +
  geom_point(data = dataset_summary %>% filter(SPLIT == TRUE),
             mapping = aes(x = RT, y = RF, text = Compound), 
             color = "#FBC438",
             size = 3, 
             alpha = 1) +
  xlim(0,25) +
  scale_y_log10(limit = c(5e17, 5e21)) +
  scale_shape(solid = FALSE) +
  labs(y = "Response factor", x = "Retention time (min)") +
  theme_light +
  annotation_logticks(sides = "l", colour = "gray")

#end----
#------using closest compound's RF----

cal_Comp_closeRT <- c()
for(i in 1:342){
  #print(dataset2[i,]$RT)
  cal_Comp_closeRT <- c(cal_Comp_closeRT, assigning_closest_cal_comp(dataset_summary[i,]$RT, cal_compounds))
}

cal_RF_closeRT <- c()
for(i in 1:342){
  #print(dataset2[i,]$RT)
  cal_RF_closeRT <- c(cal_RF_closeRT, assigning_closest_cal_comp_RF(dataset_summary[i,]$RT, cal_compounds))
}

dataset_summary <- data.frame(dataset_summary, cal_Comp_closeRT, cal_RF_closeRT)

#------using parent compound--------

parent <- parent %>%
  rename(TP = Compound,
         Compound = Parent) %>%
  left_join(dataset_summary) %>%
  rename(Parent = Compound,
         Compound = TP,
         logIE_pred_parent = logIE_pred,
         logRF_pred_parent = logRF_pred,
         RF_parent = RF) %>%
  na.omit() %>%
  select(Compound, Parent, logIE_pred_parent, logRF_pred_parent, RF_parent)

dataset_summary <- dataset_summary %>%
  left_join(parent)

#------RF dif for parents & TPs-----

dataset_summary <- dataset_summary %>%
  mutate(RF_dif_parent = case_when(
                              RF > RF_parent ~ RF/RF_parent,
                              TRUE ~ RF_parent/RF),
         RF_RTcal_comp = case_when(
                              RF > cal_RF_closeRT ~ RF/cal_RF_closeRT,
                              TRUE ~ cal_RF_closeRT/RF))

ggplot(data = dataset_summary) +
  geom_point(mapping = aes(x = RF, 
                           y = RF_parent,
                           text = Compound), 
             color = "#FBC438",
             size = 3, 
             alpha = 3/4) +
  scale_x_log10(limit = c(5e17, 5e21)) +
  scale_y_log10(limit = c(5e17, 5e21)) +
  labs(y = "Response factor parent", x = "Response factor TP") +
  geom_abline(slope = 1, intercept = 0) + 
  theme_light + 
  annotation_logticks(colour = "gray")


#---------predicting the concentrations---------
dataset <- dataset %>%
  select(-RF) %>%
  left_join(dataset_summary) %>%
  mutate(c_pred_IE = Signal_corr/(10^logRF_pred),
         c_pred_RT = Signal_corr/cal_RF_closeRT,
         c_pred_parent = Signal_corr/RF_parent) %>%
  mutate(c_pred_acc_IE =   case_when (
                            c_pred_IE > c_M ~ c_pred_IE/c_M,
                            c_pred_IE < c_M ~ c_M/c_pred_IE),
    c_pred_acc_RT =   case_when (
                            c_pred_RT > c_M ~ c_pred_RT/c_M,
                            c_pred_RT < c_M ~ c_M/c_pred_RT),
    c_pred_acc_parent =   case_when (
                            c_pred_parent > c_M ~ c_pred_parent/c_M,
                            c_pred_parent < c_M ~ c_M/c_pred_parent)) %>%
  select(Compound, SMILES, c_pred_IE, c_pred_RT, c_pred_parent, c_M, c_pred_acc_IE, c_pred_acc_RT, c_pred_acc_parent, logRF_pred, logIE_pred, RF, everything())

#-------visualizing c standards --------

#plot the predicted concentration against actual concentration
ggplot(data = dataset) +
  geom_point(mapping = aes(x =c_M, y = c_pred_IE, text = Compound), color = "#21313B", size = 3, alpha = 1/10) +
  scale_x_log10(limit = c(1e-13, 1e-8)) +
  scale_y_log10(limit = c(1e-13, 1e-8)) +
  geom_abline(slope = 1, intercept = 0, color = "#21313B",, size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Spiked concentration (M)", y = "Predicted concentration (M)") +
  theme_light  +
  annotation_logticks(colour = "gray")

ggplot(data = dataset) +
  geom_point(mapping = aes(x =c_M, y = c_pred_RT, text = Compound), color = "#21313B", size = 3, alpha = 1/10) +
  scale_x_log10(limit = c(1e-13, 1e-8)) +
  scale_y_log10(limit = c(1e-13, 1e-8)) +
  geom_abline(slope = 1, intercept = 0, color = "#21313B",, size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Spiked concentration (M)", y = "Predicted concentration (M)") +
  theme_light  +
  annotation_logticks(colour = "gray")


ggplot(data = dataset) +
  geom_point(mapping = aes(x =c_M, y = c_pred_parent, text = Compound), color = "#21313B", size = 3, alpha = 1/10) +
  scale_x_log10(limit = c(1e-13, 1e-8)) +
  scale_y_log10(limit = c(1e-13, 1e-8)) +
  geom_abline(slope = 1, intercept = 0, color = "#21313B",, size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Spiked concentration (M)", y = "Predicted concentration (M)") +
  theme_light  +
  annotation_logticks(colour = "gray")

write_delim(dataset, "Results_standards.csv", delim = ",")

#-------REAL SAMPLES------
#--------read the data --------
samples <- read_delim('Dataset_micropollutants.csv',
                       delim = ",",
                       col_names = TRUE,
                       trim_ws = TRUE)
samples <- samples %>%
  filter(samples$Sample_type == "Unknown") 

samples <- samples %>%
  left_join(dataset_summary) %>%
  left_join(Padel_data)

samples <- samples %>%
  mutate(Spiked = case_when(
           str_detect(samples$Sample_name, "spike") ~ "yes",
           str_detect(samples$Sample_name, "STD") ~ "yes",
           TRUE ~ "no")) %>%
  select(Compound, Signal, ISTD_response, ISTD_area_mean, Sample_type, everything()) %>%
  mutate(Signal_corr = Signal/(ISTD_response/ISTD_area_mean)) %>%
  select(Compound, Signal_corr, everything())

#--------calculating concentration-------
samples <- samples %>%
  mutate(
    c_M = Concentration_sample* (10^-9)/MW, #the starting concentration is in ng/l
    Signal_corr = Signal_corr*IC
  ) %>%
  filter(c_M > 0) %>%
  mutate(c_pred_IE = Signal_corr/(10^logRF_pred),
         c_pred_RT = Signal_corr/cal_RF_closeRT,
         c_pred_parent = Signal_corr/RF_parent) %>%
  mutate(c_pred_acc_IE =   case_when (
                                    c_pred_IE > c_M ~ c_pred_IE/c_M,
                                    c_pred_IE < c_M ~ c_M/c_pred_IE),
    c_pred_acc_RT =   case_when (
                                    c_pred_RT > c_M ~ c_pred_RT/c_M,
                                    c_pred_RT < c_M ~ c_M/c_pred_RT),
    c_pred_acc_parent =   case_when (
                                    c_pred_parent > c_M ~ c_pred_parent/c_M,
                                    c_pred_parent < c_M ~ c_M/c_pred_parent)) %>%
  select(Compound, SMILES, c_pred_IE, c_pred_RT,c_pred_parent,  c_M, c_pred_acc_IE, c_pred_acc_RT,c_pred_acc_parent, logRF_pred, logIE_pred, RF, everything())

ggplot(data = samples%>% filter(samples$Spiked == "no")) +
  geom_point(mapping = aes(x =c_M, y = c_pred_RT, text = Compound), color = "#21313B", size = 3, alpha = 2/3) +
  geom_point(mapping = aes(x =c_M, y = c_pred_parent, text = Compound), color = "#ED3A53", size = 3, alpha = 2/3) +
  geom_point(mapping = aes(x =c_M, y = c_pred_IE, text = Compound), color = "#FBC438", size = 3, alpha = 2/3) +
  facet_wrap(~ Sample_name, ncol = 8) +
  scale_x_log10(limit = c(1e-13, 1e-8)) +
  scale_y_log10(limit = c(1e-13, 1e-8)) +
  geom_abline(slope = 1, intercept = 0, color='steelblue', size = 1) +
  scale_shape(solid = FALSE) +
  labs(x = "Measured concentration (M)", y = "Predicted concentration (M)") +
  theme_light_wide +
  theme(axis.text.x = element_text(angle = -90, hjust = 1)) +
  annotation_logticks(colour = "gray")

write_delim(samples, "Results_samples.csv", delim = ",")