# Loading in libraries ----------------------------------------------------
#Data processing libraries
library(tidyverse)
library(data.table)
library(DescTools)
library(readxl)

# statistical libraries
library(randomForest)
library(lme4)
library(MuMIn)

#Setting specific functions that overlap with others to guaruntee correct imputation
map <- purrr::map
select <- dplyr::select
tidy <- broom::tidy
rename <- dplyr::rename
mutate <- dplyr::mutate

## Defining some functions that I use for data visualization
zscore <- function(x) {
  (x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

gen_theme <-  function(x){
  theme(plot.margin = unit(c(1,1,1.5,1.2), 'cm'),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 25),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
        panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"))
}


# Loading in raw dataframes -----------------------------------------------
# Reading in all of the raw data for metabolomics, genotypes, and proteomics
# Then a new column is created to designate the type of data so that we can make each dataframe a tibble and bind them together
# rawMetab <- read_csv('~/Documents/GitHub/yuHCT/data/raw/metabolomics_20240415.csv')%>%
#   mutate(dataType = 'metabolomics')

metabIIDs <- read_xlsx('~/Documents/GitHub/yuHCT/data/raw/ID_Key_iid_and_CDP_ID.xlsx')%>%
  select(8:9)%>%
  rename(sampleName = 1, 
         iid = 2)%>%
  mutate(sampleName = gsub('CDP_', '', sampleName))

badNames <- c('80_OHT_AM',
              '83_PSM_AM',
              '84_MJE_AM',
              '86_ECC_AM',
              '87_OTR_AM')


rawMetab <- read_xlsx('~/Documents/GitHub/yuHCT/data/raw/2019.04.22.CDP_Final_Data_Extraction.xlsx', sheet =2)%>%
  select(9:ncol(.))%>%
  rownames_to_column(var = 'metabolite')%>%
  mutate(m = 'M')%>%
  unite(metabolite, c('m', 'metabolite'), sep = "")%>%
  pivot_longer(2:ncol(.), names_to = 'sampleName', values_to = 'xic')%>%
  mutate(xic = log10(xic + 1))%>%
  pivot_wider(names_from = 'metabolite', values_from = 'xic')%>%
  filter(!sampleName %in% badNames)%>%
  left_join(metabIIDs, by = 'sampleName')%>%
  # select(-sampleName)%>%
  mutate(dataType = 'metabolomics')%>%
  select(iid, dataType, everything())

metadataMetab <- read_csv('~/Documents/GitHub/yuHCT/data/raw/103019_metabolite_explain.csv')

rawGenotype <- read_csv('~/Documents/GitHub/yuHCT/data/raw/genotype_20240415.csv')%>%
  mutate(dataType = 'genotype')
# group_by(dataType)%>%
# nest()

rawPhenotype <- read_csv('~/Documents/GitHub/yuHCT/data/raw/phenotype_20240603.csv')%>%
  select(-X23)

rawProteomics <- read_csv('~/Documents/GitHub/yuHCT/data/raw/proteomics_20240415.csv')%>%
  mutate(dataType = 'proteomics')
# group_by(dataType)%>%
# nest()


# Here we have forced each data type into tibbles and put each tibble as a single row.
# This allows us to map the same model (code) over each data set separately
joinAndNest <- function(x) {
  left_join(x, rawPhenotype%>% 
              select(iid, hct:hchvr, age, sex, bmi, smoker),
            by = 'iid')%>%
    select(hct:hchvr, iid, everything())%>%
    pivot_longer(1:17, names_to = 'responseVariable', values_to = 'responseValue')%>%
    group_by(dataType, responseVariable)%>%
    nest()%>%
    mutate(data = map(data, ~ select(.x, iid, age, sex, bmi, smoker, responseValue, everything())%>%
                        filter(!is.na(responseValue))))
}

## Applying the join and Nest function to make sure they are all run the same way
nestedProteomics <- rawProteomics%>%
  joinAndNest()%>%
  mutate(data = map(data, ~pivot_longer(.x, 8:ncol(.), names_to = 'feature', values_to = 'values')))

nestedGenoType <- rawGenotype%>%
  joinAndNest()%>%
  mutate(data = map(data, ~pivot_longer(.x, 8:ncol(.), names_to = 'feature', values_to = 'values')%>%
                      mutate(values = as.factor(values))))


nestedMetabolomes <- rawMetab%>%
  select(-sampleName)%>%
  joinAndNest()

# Finding a gap fill value to replace NA's in metabolomics data. 
# Value = min XIC / 2
gapFillValue <- (rawMetab%>% 
                   select(-sampleName)%>%
                   select(iid, dataType, everything())%>%
                   pivot_longer(3:ncol(.), names_to = 'metabolite', values_to = 'xic')%>%
                   select(xic)%>%
                   filter(!is.na(xic))%>%
                   as.vector()%>%
                   min())/2


# Random forest values
## These were originally derived in the below code but to save time and computational power were saved and are reuploaded each run

setwd('~/Documents/GitHub/yuHCT/data/analysis/randomForest/')

responses <- dir(path = '~/Documents/GitHub/yuHCT/data/analysis/randomForest/', pattern = '.csv')

rfValues <- responses%>%
  map(read_csv)%>%
  tibble()%>%
  mutate(dataType = 'metabolomics',
         responseVariable = responses)%>%
  rename(data = 1)%>%
  mutate(responseVariable = gsub('.csv', '', responseVariable))

setwd('~/Documents/GitHub/yuHCT/data/')
  
# Cleaning -- randomForest dataframes ---------------------------------------------------
combinedData <- nestedMetabolomes%>%
  mutate(rForest = map(data, ~select(.x, -c(iid, age, sex, bmi, smoker))%>%
                         mutate_at(2:ncol(.), ~ifelse(is.na(.), gapFillValue,.)))) # This imputes the gap fill values for any NA's. We probably want to actually just filter out rows with NA's

# STATS -- set seed -------------------------------------------------------
set.seed(12683018)


# STATS -- qq plots ------------------------------------------------------
pdf('~/Documents/GitHub/yuHCT/data/plots/qqPlots.pdf', width = 12, height = 10)
rawPhenotype%>%
  select(iid, hct:hchvr)%>%
  pivot_longer(2:ncol(.), names_to = 'responseVariable', values_to = 'responseValue')%>%
  group_by(responseVariable)%>%
  nest()%>%
  mutate(data = map2(data, responseVariable, ~ pull(.x, responseValue)%>%
                      car::qqPlot() +
                       title(.y)))
dev.off()
# STATS -- random Forest --------------------------------------------------
#Random forest 

## This section is how the random forest was originally computed for all of the phenotypes. However, because it takes so long to run
## the random forest, we now just import the analysis files belows.
# # rfDataFrame <- combinedData%>%
# #   mutate(rForest = map(rForest, ~ randomForest(responseValue ~ ., .x,importance = TRUE, proximity = TRUE, ntree = 100000)))
# 
# rfValues <- rfDataFrame%>%
#   mutate(rForest = map(rForest, ~ (.x$importance%>%
#                                      as.data.frame()%>%
#                                      rownames_to_column("feature")%>%
#                                      # mutate(feature = gsub("X", "", feature),
#                                      #        feature = gsub("\\.", "-", feature))%>%
#                                      filter(`%IncMSE` >= mean(`%IncMSE`) + 2*sd(`%IncMSE`)))$feature%>%
#                          as.vector()))%>%
#     mutate(data = map2(data, rForest, ~pivot_longer(.x, 8:ncol(.), names_to = 'feature', values_to = 'values')%>%
#                          filter(feature %in% .y)),
#            name = str_c('~/Documents/GitHub/yuHCT/data/analysis/randomForest', responseVariable, '.csv'))%>%
#     select(-rForest)

# map2(rfValues$data, rfValues$name, ~write_csv(.x%>% as.data.frame(), .y))


rfValues%>%
  # select(-data)%>%
  unnest(data)%>%
  pull(feature)%>%
  unique()%>%
  length()


# STATS -- lmer -----------------------------------------------------------
lmerData <- rfValues%>%
  bind_rows(nestedProteomics,
            nestedGenoType)%>%
  mutate(lmer = map(data, ~group_by(.x,feature)%>%
                      filter(!feature %in% c('M1134', 'M1662', 'M3801', 
                                             'M2055', 'M2465', 'M3453', 
                                             'M3678', 'M3759', 'M4054', 
                                             'M4145', 'M5009', 'M5370',
                                             'M5479', 'M5604', 'M5739',
                                             'M5905', 'M6663', 'M7389',
                                             'M788', 'M8052', 'M8432', 'M999')  
                      )%>% 
                      mutate(ageRange = case_when(age < 30 ~ '18-29',
                                             age >= 30 & age < 40 ~ '30-439',
                                             age >= 40 & age < 50 ~ '40-49',
                                             age >= 50 & age < 60 ~ '50-59',
                                             age >=60 & age < 70 ~ '60-69'),
                             bmiRange = case_when(bmi < 20 ~ '<20',
                                                  bmi >= 20 & bmi < 25 ~'20-25',
                                                  bmi >= 25 & bmi < 30 ~'25-30',
                                                  bmi >= 30 & bmi < 35 ~'30-35',
                                                  bmi >= 35 ~ '>35'))%>%# All of these metabolites had too few observations, or grouping variables (sex/smoker/etc) were singular
                      nest()%>%
                      mutate(data = map(data, ~ filter(.x, !is.na(values),
                                                       !is.na(sex),
                                                       !is.na(smoker),
                                                       !is.na(ageRange),
                                                       !is.na(bmiRange))%>%
                                          lmer(responseValue~values + (1|sex) + (1|smoker) + (1|ageRange) + (1|bmiRange), ## All of these controls limit the n ~ 55-70
                                               data =., control =  lmerControl(check.nlev.gtr.1 = "ignore",
                                                                               check.conv.singular = 'ignore',
                                                                               check.nobs.vs.nRE = 'ignore'))),
                             lmerPvalues = map(data, ~ car::Anova(.x)%>%
                                                 rownames_to_column(var = 'var')%>%
                                                 select(var, `Pr(>Chisq)`)),
                             lmerR2 = map(data, ~ r.squaredGLMM(.x)%>%
                                            tidy()))))%>%
  select(-c(data))%>%
  unnest(lmer)%>%
  unnest(c(lmerPvalues, lmerR2))%>%
  mutate(FDR = p.adjust(`Pr(>Chisq)`, method ='BH'))%>%
  filter(FDR < 0.05)


significantLmer <- lmerData%>%
  select(-data)

write_csv(significantLmer, '~/Documents/GitHub/yuHCT/data/analysis/allPvalues.csv')

significantMetabolites <- significantLmer%>% 
  filter(dataType == 'metabolomics')%>% 
  group_by(responseVariable)%>%
  mutate(min= min(R2m), 
         max = max(R2m), 
         maxP = max(FDR), 
         n = 1, 
         n = sum(n))%>% 
  summarize_if(is.numeric, mean)%>% 
  select(responseVariable, min, max, maxP, n)

write_csv(significantMetabolites, '~/Documents/GitHub/yuHCT/data/analysis/significantMetaboliteSummary.csv')

# Hierarchical clustering -------------------------------------------------
# circosFeatureClusters <- significantLmer%>%
#   ungroup()%>%
#   filter(dataType == 'metabolomics',
#          R2m > 0.1)%>%
#   select(responseVariable, feature)%>%
#   unique()%>%
#   group_by(feature)%>%
#   summarise(cluster = toString(responseVariable))%>%
#   mutate(cluster = gsub(', ', '_', cluster))
# 
# 
# hcFeatures <- (significantLmer%>%
#   ungroup()%>%
#   filter(dataType == 'metabolomics',
#          R2m > 0.1))$feature%>%
#   as.vector()%>%
#   unique()
# 
# 
# hcMetab <- rawMetab%>%
#   filter(sampleName != 'V120')%>%
#   select(-c(sampleName, dataType))%>% ##sampleName needs to be changed to iid when James gives the correct format
#   pivot_longer(2:ncol(.), names_to = 'features', values_to = 'xic')%>%
#   filter(features %in% hcFeatures)%>%
#   group_by(features)%>%
#   mutate(xic = ifelse(is.na(xic), gapFillValue, xic),
#          xic = zscore(xic))%>%
#   ungroup()%>%
#   pivot_wider(names_from = 'features', values_from = 'xic')%>%
#   left_join(rawPhenotype%>%
#               select(-c(age:smoker))%>%
#               pivot_longer(2:ncol(.), names_to = 'phenotype', values_to = 'value')%>%
#               group_by(phenotype)%>%
#               mutate(value = zscore(value))%>%
#               pivot_wider(names_from = 'phenotype', values_from = 'value'),
#             by = 'iid')%>%
#   column_to_rownames(var = 'iid')%>%
#   pheatmap::pheatmap(color = brewer.pal(n = 9, name = "Greys"), 
#                      fontsize_col = 3, fontsize_row = 4, cutree_cols = 14)
# 
# #hc Orders
# hcOrders <- tibble(feature = hcMetab$tree_col[['labels']], 
#                    order =  hcMetab$tree_col[['order']])%>%
#   arrange(order)
# 
# clusters <- cutree(hcMetab$tree_col, k = 14)%>%
#   as.data.frame()%>%
#   rownames_to_column(var = 'feature')%>%
#   rename(cluster = 2)
# 
# 
# #HC plots
# pdf('~/Documents/GitHub/yuHCT/data/plots/hcMetabolites.svg', width = 15, height = 10)
# hcMetab
# dev.off()


# hc -- phenotype hc ------------------------------------------------------
## Once I have the new IID Dataframe I can change sampleName out for it and then this should work a lot smoother as a subplot of the HC
# hcIidOrder <- (tibble(feature = hcMetab$tree_row[['labels']], 
#                      order =  hcMetab$tree_row[['order']])%>%
#   arrange(order))$feature
# 
# test <- hcMetab$tree_row$order
# 
# hcPhenotype <- rawPhenotype%>% 
#   select(iid, hct, hb, SpO2, sbp, dbp, glucose, insulin, iron, hvr)%>%
#   pivot_longer(2:ncol(.), names_to = 'phenotype', values_to = 'values')%>%
#   left_join(metabIIDs, by = 'iid')%>%
#   filter(!sampleName %in% badNames)%>%
#   group_by(phenotype)%>%
#   mutate(iid = as.factor(iid),
#          iid = fct_relevel(iid, hcIidOrder),
#          zscore = zscore(values))%>%
#   ungroup()%>%
#   select(-c(sampleName, values))%>%
#   # ggplot() +
#   # geom_tile(aes(phenotype, iid, fill = zscore))
#   pivot_wider(names_from = 'phenotype', values_from = 'zscore')%>%
#   column_to_rownames(var = 'iid')%>%
#   pheatmap::pheatmap(cluster_rows = FALSE, na.omit = TRUE, fontsize_row = 4, display_numbers = test,
#                      color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
# dev.off()  
# 
# 
# # hcPhenotype$tree_row <- hcMetab$tree_row
# 
# 
# 
# 
# pdf('~/Documents/GitHub/yuHCT/data/plots/hcPhenotypes.pdf', width = 15, height = 10)
# hcPhenotype
# dev.off()


# Circos Plot -- Chromosomes and bands -------------------------------------------------
# Phenotypes have to be created separately and bound to the others
phenotypeChromosome <- significantLmer%>%
  filter(R2m > 0.1)%>%
  ungroup()%>%
  select(responseVariable)%>%
  unique()%>%
  mutate(chr = 'chr -',
         dataType = 'Phenotype',
         color = 'vlgrey',
         zero = 0,
         end = 30,
         dataType2 = dataType)%>%
  group_by(chr, dataType, dataType2, color)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  unite(column, c('chr', 'dataType', 'dataType2', 'zero', 'end', 'color'), sep = ' ')

#Phenotype bands
phenotypeBands <- significantLmer%>%
  filter(R2m > 0.1)%>%
  ungroup()%>%
  select(responseVariable)%>%
  unique()%>%
  mutate(chr = 'band',
         dataType = 'Phenotype',
         color = 'vlgrey',
         zero = 0,
         end = 30,
         dataType2 = dataType)%>%
  rename(feature = responseVariable)%>%
  mutate(bandEnd = cumsum(end),
         num = 1,
         num = cumsum(num),
         start = ifelse(num>1, lag(bandEnd), 0),
         bandEnd = format(as.numeric(bandEnd), scientific = FALSE, trim = TRUE),
         start = format(as.numeric(start), scientific = FALSE, trim = TRUE))%>%
  mutate(feature2 = feature)%>%
  select(-c(zero, num, end, dataType2))%>%
  unite(column, c('chr', 'dataType', 'feature', 'feature2', 'start', 'bandEnd', 'color'), sep = ' ')
  

## Main prep for the rest of the circos Diagram
circosPrep <- significantLmer%>%
  filter(R2m > 0.1)%>%
  mutate(chr = 'chr -',
         zero = 0,
         # end = 100,
         end = ifelse(dataType != 'metabolomics', 30, 10),
         feature2 = feature)%>%
  select(chr, dataType, feature, feature2, zero, end, responseVariable, R2m)%>%
  mutate(dataType = case_when(dataType == 'phenotype' ~ 'Phenotype',
                              dataType == 'metabolomics' ~ 'Metabalome',
                              dataType == 'genotype' ~ 'Genotype',
                              dataType == 'proteomics' ~ 'Proteome'))

#Creating the four chromosomes
circosChromosome <- circosPrep%>%
  select(-c(R2m))%>%
  mutate(color = case_when(dataType == 'Phenotype' ~ 'vlgrey',
                           dataType == 'Metabalome' ~ 'vlgreen',
                           dataType == 'Genotype' ~ 'vlblue',
                           dataType == 'Proteome' ~ 'vlred'))%>%
  ungroup()%>%
  select(chr, dataType, feature, zero, color, zero, end, color)%>%
  unique()%>%
  group_by(chr, dataType, zero, color)%>%
  summarize_if(is.numeric, sum)%>%
  mutate(dataType2 = dataType)%>%
  unite(column, c('chr', 'dataType', 'dataType2', 'zero', 'end', 'color'), sep = ' ')%>%
  bind_rows(phenotypeChromosome)

# Creating the badnds
circosBands <- circosPrep%>%
  mutate(chr = 'band')%>%
  ungroup()%>%
  select(chr, dataType, feature, zero, end)%>%
  unique()%>%
  group_by(dataType)%>%
  left_join(circosFeatureClusters, by = 'feature')%>%
  arrange(cluster)%>%
  mutate(bandEnd = cumsum(end),
         num = 1,
         num = cumsum(num),
         start = ifelse(num>1, lag(bandEnd), 0),
         color = case_when(dataType == 'Phenotype' ~ 'vlgrey',
                           dataType == 'Metabalome' ~ 'vlgreen',
                           dataType == 'Genotype' ~ 'vlblue',
                           dataType == 'Proteome' ~ 'vlred'),
         bandEnd = format(as.numeric(bandEnd), scientific = FALSE, trim = TRUE),
         start = format(as.numeric(start), scientific = FALSE, trim = TRUE))%>%
  mutate(feature2 = feature)%>%
  select(-c(zero, num, end, cluster))%>%
  unite(column, c('chr', 'dataType', 'feature', 'feature2', 'start', 'bandEnd', 'color'), sep = ' ')%>%
  bind_rows(phenotypeBands)

#Combining into a single csv to write as txt
circosBinaryPrint <- circosChromosome%>%
  bind_rows(circosBands)%>%
              # select(-c(responseVariable, lmerR2)))%>%
  as.matrix()

write_lines(circosBinaryPrint, '~/Documents/GitHub/yuHCT/data/plots/circos/circos.txt')



# Circos Plot -- Tiles ----------------------------------------------------
circosTiles <- circosBands%>%
  separate(column, c('chr', 'name', 'network', 'network2', 'start', 'end', 'color'), sep = ' ')%>%
  # left_join(clusters, by = c('network'='feature'))%>%
  mutate(colorVal = 'color=',
         color = case_when(color=='vlgrey'~ 'lgrey',
                           color=='vlgreen'~ 'lgreen',
                           color=='vlred'~ 'lred',
                           color=='vlblue'~ 'lblue'))%>%
  select(name, start, end, color, colorVal)%>%
  unite(color, c(colorVal, color), sep = '')%>%
  unite(column, c(name, start, end, color), sep = ' ')%>%
  as.matrix()

write_lines(circosTiles, '~/Documents/GitHub/yuHCT/data/plots/circos/tiles.txt')
  

# Circos Plot -- Links ----------------------------------------------------
variableLinks <- circosPrep%>%
  ungroup()%>%
  select(feature, responseVariable, R2m)%>%
  left_join(circosBands%>%
              separate(column, c('chr', 'name', 'feature', 'feature2', 'start', 'bandEnd', 'color'), sep = ' ')%>%
              select(name, feature, start, bandEnd)%>%
              unique(), 
            by = 'feature')%>%
  ungroup()%>%
  group_by(feature)%>%
  mutate(bandStart = as.numeric(start),
         bandEnd = as.numeric(bandEnd),
         n = 1,
         sumN = sum(n),
         num = cumsum(n),
         size = ifelse(name == 'Metabalome', 6, 26),
         bandsize = floor(size/sumN),
         bandsize = ifelse(num == 1 & num != sumN, 0, bandsize),
         linkStart = ifelse(num > 1, cumsum(bandsize) + bandStart + 2, bandStart + 2),
         linkEnding = ifelse(num < sumN, lead(linkStart), linkStart+bandsize),
         colorBin = 'color=',
         color = case_when(#Metabalome
           R2m <= 0.2 & name == 'Metabalome' ~ 'vvlgreen',
           R2m > 0.2 & R2m <= 0.3 & name == 'Metabalome' ~ 'vlgreen',
           R2m > 0.3 & R2m <= 0.4 & name == 'Metabalome' ~ 'lgreen',
           R2m > 0.4 & R2m <= 0.5 & name == 'Metabalome' ~ 'green',
           R2m > 0.5 & R2m <= 0.6 & name == 'Metabalome' ~ 'dgreen',
           R2m > 0.6 & R2m <= 0.7 & name == 'Metabalome' ~ 'vdgreen',
           R2m > 0.7 & R2m <= 0.8 & name == 'Metabalome' ~ 'vvdgreen',
           #Proteome colors
           R2m <= 0.13 & name == 'Proteome' ~ 'vvlred',
           R2m > 0.13 & R2m <= 0.16 &  name == 'Proteome' ~ 'vlred',
           R2m > 0.16 & R2m <= 0.19 & name == 'Proteome' ~ 'lred',
           R2m > 0.19 & R2m <= 0.22 & name == 'Proteome' ~ 'red',
           R2m > 0.22 & R2m <= 0.25 & name == 'Proteome' ~ 'dred',
           R2m > 0.25 & R2m <= 0.28 & name == 'Proteome' ~ 'vdred',
           R2m > 0.28 & R2m <= 0.33 & name == 'Proteome' ~ 'vvdred',
           #Genotype colors
          R2m <= 0.15 & name == 'Genotype' ~ 'vvlblue',
           R2m > 0.15 & R2m <= 0.18 & name == 'Genotype' ~ 'vlblue',
           R2m > 0.18 & R2m <= 0.21 & name == 'Genotype' ~ 'lblue',
           R2m > 0.21 & R2m <= 0.24 & name == 'Genotype' ~ 'blue',
           R2m > 0.24 & R2m <= 0.27 & name == 'Genotype' ~ 'dblue',
           R2m > 0.27 & R2m <= 0.30 & name == 'Genotype' ~ 'vdblue',
           R2m > 0.30 & R2m <= 0.35 & name == 'Genotype' ~ 'vvdblue'
         ))%>%
  ungroup()%>%
  unite(color, c('colorBin', 'color'), sep = '')%>%
  rownames_to_column(var = 'linkID')%>%
  unite(column, c(linkID, name, linkStart, linkEnding, color), sep = ' ')%>% 
  select(column)
  

responseLinks <- circosPrep%>%
  ungroup()%>%
  select(feature, dataType, responseVariable, R2m)%>%
  left_join(phenotypeBands%>%
              separate(column, c('chr', 'name', 'responseVariable', 'feature2', 'start', 'bandEnd', 'color'), sep = ' ')%>%
              select(name, responseVariable, start, bandEnd), 
            by = 'responseVariable')%>%
  ungroup()%>%
  group_by(responseVariable)%>%
  mutate(bandStart = as.numeric(start),
         bandEnd = as.numeric(bandEnd),
         n = 1,
         sumN = sum(n),
         num = cumsum(n),
         bandsize = floor(26/sumN),
         bandsize = ifelse(num == 1 & num != sumN, 0, bandsize),
         linkStart = ifelse(num > 1, cumsum(bandsize) + bandStart + 2, bandStart + 2),
         linkEnding = ifelse(num < sumN, lead(linkStart), linkStart+bandsize),
         colorBin = 'color=',
         color = case_when(#Metabalome
           R2m <= 0.2 & dataType == 'Metabalome' ~ 'vvlgreen',
           R2m > 0.2 & R2m <= 0.3 & dataType == 'Metabalome' ~ 'vlgreen',
           R2m > 0.3 & R2m <= 0.4 & dataType == 'Metabalome' ~ 'lgreen',
           R2m > 0.4 & R2m <= 0.5 & dataType == 'Metabalome' ~ 'green',
           R2m > 0.5 & R2m <= 0.6 & dataType == 'Metabalome' ~ 'dgreen',
           R2m > 0.6 & R2m <= 0.7 & dataType == 'Metabalome' ~ 'vdgreen',
           R2m > 0.7 & R2m <= 0.8 & dataType == 'Metabalome' ~ 'vvdgreen',
           #Proteome colors
           R2m <= 0.13 & dataType == 'Proteome' ~ 'vvlred',
           R2m > 0.13 & R2m <= 0.16 &  dataType == 'Proteome' ~ 'vlred',
           R2m > 0.16 & R2m <= 0.19 & dataType == 'Proteome' ~ 'lred',
           R2m > 0.19 & R2m <= 0.22 & dataType == 'Proteome' ~ 'red',
           R2m > 0.22 & R2m <= 0.25 & dataType == 'Proteome' ~ 'dred',
           R2m > 0.25 & R2m <= 0.28 & dataType == 'Proteome' ~ 'vdred',
           R2m > 0.28 & R2m <= 0.33 & dataType == 'Proteome' ~ 'vvdred',
           #Genotype colors
           R2m <= 0.15 & dataType == 'Genotype' ~ 'vvlblue',
           R2m > 0.15 & R2m <= 0.18 & dataType == 'Genotype' ~ 'vlblue',
           R2m > 0.18 & R2m <= 0.21 & dataType == 'Genotype' ~ 'lblue',
           R2m > 0.21 & R2m <= 0.24 & dataType == 'Genotype' ~ 'blue',
           R2m > 0.24 & R2m <= 0.27 & dataType == 'Genotype' ~ 'dblue',
           R2m > 0.27 & R2m <= 0.30 & dataType == 'Genotype' ~ 'vdblue',
           R2m > 0.30 & R2m <= 0.35 & dataType == 'Genotype' ~ 'vvdblue'
           ))%>%
  ungroup()%>%
  unite(color, c('colorBin', 'color'), sep = '')%>%
  rownames_to_column(var = 'linkID')%>%
  unite(column, c(linkID, name, linkStart, linkEnding, color), sep = ' ')%>% 
  select(column)
  

linkJoiner <- variableLinks%>%
  bind_rows(responseLinks)%>%
  as.matrix()

write_lines(linkJoiner, '~/Documents/GitHub/yuHCT/data/plots/circos/circoslinks.txt')



# CircosPlot -- text labels -----------------------------------------------
text <- circosBands%>%
  separate(column, c('chr', 'name', 'network', 'network2', 'start', 'end', 'color'), sep = ' ')%>%
  select(name, start, end, network, color)%>%
  mutate(color = 'color=vvdgrey')%>%
  unite(column, c(name, start, end, network, color), sep = ' ')%>%
  as.matrix()

write_lines(text, '~/Documents/GitHub/yuHCT/data/plots/circos/text.txt')

# Cluster text for building cluster mapping
clusterText <- circosBands%>%
  separate(column, c('chr', 'name', 'feature', 'network2', 'start', 'end', 'color'), sep = ' ')%>%
  left_join(circosFeatureClusters, by = 'feature')%>%
  select(name, start, end, cluster, color)%>%
  mutate(color = 'color=vvdgrey')%>%
  unite(column, c(name, start, end, cluster, color), sep = ' ')%>%
  as.matrix()

write_lines(clusterText, '~/Documents/GitHub/yuHCT/data/plots/circos/clusterText.txt')



# CircosPLot -- Lines for each cluster ------------------------------------
clusterBands <- circosBands%>%
  separate(column, c('chr', 'name', 'feature', 'network2', 'start', 'end', 'color'), sep = ' ')%>%
  left_join(circosFeatureClusters, by = 'feature')%>%
  mutate(val = 1)%>%
  select(name, start, end, val, color)%>%
  mutate(color = 'color=vvdgrey')%>%
  unite(column, c(name, start, end, val, color), sep = ' ')%>%
  as.matrix()

write_lines(clusterBands, '~/Documents/GitHub/yuHCT/data/plots/circos/lines.txt')
  








