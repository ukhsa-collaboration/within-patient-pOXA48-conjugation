library( MASS)                        # for glm.nb()
library(ggplot2)
library(tidyverse)

my_data <- read.csv("data/genetic_distance_table.csv")

# Remove the comparison to reference plasmids
my_data <- my_data[my_data$pl1_pat != "Reference", ]
my_data <- my_data[my_data$pl2_pat != "Reference", ]

# create a column which checks if the plasmids is being compared to itself
my_data$same_sample = my_data$pl1_masked == my_data$pl2_masked

# convert within patient into a 0 for not within the same patient and 1 for within the same patient
my_data$population <- as.numeric(my_data$within_patient == "True")

# Filter out comparisons to the same sample
modelled_df <- my_data[my_data$same_sample ==FALSE, c("SNPs", "population")]

# Fit negative binomial model
glmFitNB <- glm.nb(SNPs ~ population, data=modelled_df)
summary(glmFitNB)                     # negative binomial model

# extract coefficients
est <- cbind(Estimate = coef(glmFitNB), confint(glmFitNB))

est

#  incident rate ratios
exp(est)


## prediction

newdata1 <- data.frame(population = c(0, 1))
newdata1$snp <- predict(glmFitNB, newdata1, type = "response")
newdata1

##
newdata2 <- data.frame(
  population = c(0, 1))
newdata2 <- cbind(newdata2, predict(glmFitNB, newdata2, type = "link", se.fit=TRUE))
newdata2 <- within(newdata2, {
  SNPs <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

# ggplot(newdata2, aes(population, SNPs)) +
#   geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25) +
#   geom_line( size = 2) +
#   labs(x = "POPULATION", y = "SNPS")

ggplot(newdata2, aes(population, SNPs)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = LL, ymax = UL), size = 0.5, width = 0.1) + 
  scale_x_continuous(breaks = c(0, 1))

# histograms to display data
hist(my_data[my_data$population == 1,]$SNPs, breaks = 100)
hist(my_data[my_data$population == 0,]$SNPs, breaks = 100)





#################### REPEATING ANALYSIS FOR NORTH WEST ##########################
#################################################################################

stop("You will need the unmasked data with SAMPLE_ID-ID's to run this next bit of the script")

# Read in epi data
epi_data <- read.csv("data/supplementary_data.csv")

# Upload unmasked data
unmasked_data <- read.csv("data/genetic_distance_table.csv")

# Remove the comparison to reference plasmids
unmasked_data <- unmasked_data[unmasked_data$pl1_pat != "Reference", ]
unmasked_data <- unmasked_data[unmasked_data$pl2_pat != "Reference", ]
unmasked_data <- unmasked_data[unmasked_data$pl1_pat != "JN626286", ]
unmasked_data <- unmasked_data[unmasked_data$pl2_pat != "JN626286", ]

# create a column which checks if the plasmids is being compared to itself
unmasked_data$same_sample = unmasked_data$pl1 == unmasked_data$pl2

# convert within patient into a 0 for not within the same patient and 1 for within the same patient
unmasked_data$population <- as.numeric(unmasked_data$within_patient == "True")

# add region data:
region_data_1 <- epi_data %>% select(SAMPLE_ID, REGION) %>% rename(pl1 = SAMPLE_ID, pl1_region = REGION)
region_data_2 <- epi_data %>% select(SAMPLE_ID, REGION) %>% rename(pl2 = SAMPLE_ID, pl2_region = REGION)
unmasked_data <- left_join(unmasked_data, region_data_1, by = "pl1")
unmasked_data <- left_join(unmasked_data, region_data_2, by = "pl2")

# Filter out comparisons to the same sample & Separate into NW and NOT NW
NW_modelled_df <- unmasked_data %>% filter(same_sample == FALSE) %>% filter(pl1_region == "N WEST") %>% filter(pl2_region == "N WEST") %>% select(SNPs, population)
NEGATION_NW_modelled_df <- unmasked_data %>% filter(same_sample == FALSE) %>% filter(pl1_region != "N WEST" | pl2_region != "N WEST") %>% select(SNPs, population)
NOT_NW_modelled_df <- unmasked_data %>% filter(same_sample == FALSE) %>% filter(pl1_region != "N WEST") %>% filter(pl2_region != "N WEST") %>% select(SNPs, population)

# Modelling N WEST:
# Fit negative binomial model
NW_glmFitNB <- glm.nb(SNPs ~ population, data=NW_modelled_df)
summary(NW_glmFitNB)                     # negative binomial model
# extract coefficients
NW_est <- cbind(Estimate = coef(NW_glmFitNB), confint(NW_glmFitNB))
NW_est
#  incident rate ratios
exp(NW_est)

# TO GET CONFIDENCE INTERVAL FOR DISPERSION PARAMETER: Estimate +- 1.96 * Standard Error
