devtools::install_github("Teichlab/sctkr")
library(sctkr)

install.packages("tidyverse")
library(tidyverse)

meta_for_GLMM <- read.csv("~/ALI-SARS/adata_fig_2_metadata_glmm.csv")

#separate by age group
meta_for_GLMM_paed <- meta_for_GLMM[(meta_for_GLMM$age_group=="Paediatric") , ]
meta_for_GLMM_adult <- meta_for_GLMM[(meta_for_GLMM$age_group=="Adult") , ]
meta_for_GLMM_elderly <- meta_for_GLMM[(meta_for_GLMM$age_group=="Elderly") , ]

order = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1',
          'Basaloid-like 2', 'Cycling basal', 'Hillock','Goblet 1',
          'Goblet 2', 'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+',
          'Goblet 2 inflammatory', 'Secretory', 'Secretory 2',
          'Secretory 3', 'Secretory 4', 'Transit epi', 'Transit epi 2', 'Ciliated 1',
          'Ciliated 2', 'Deutorosomal', 'Ionocyte', 'Squamous')

rev_order <- rev(order)

# all --------------------------------------------------------------
result <- CellTypeCompositionAnalysis(
  meta_for_GLMM, 
  'souporcell_cluster', 
  'annotation_v5', 
  colVarCats=c(
    'donor_id', 'treatment', 'age_group', "kit_version", "spike.in_primer", 'pool', 'time', 'age_treatment'
  ), 
  colVarNums = c(), 
  save = 'ALI.CTCA'
)

print(plot_ranef(
  result$ranef,
  vars = list( 
   age_treatment =  c(
      'Paediatric-mock', 
      'Paediatric-SARS',
      'Adult-mock', 
      'Adult-SARS',
      'Elderly-mock', 
      'Elderly-SARS')), 
   celltype_order = rev_order, 
  maxFC=2,  
) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))


print(plot_ranef(
  result$ranef,
  vars = list(treatment = c('mock', 'SARS')), 
  references = list(treatment = 'mock'),
  celltype_order = rev_order, 
  maxFC=3, 
) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + ggtitle("Paediatric"))


print(plot_sdse(result$sdse, 'souporcell_cluster'))

# peadiatric --------------------------------------------------------------
result <- CellTypeCompositionAnalysis(
  meta_for_GLMM_paed, 
  'souporcell_cluster', 
  'annotation_v5', 
  colVarCats=c(
    'donor_id', 'treatment', "kit_version", "spike.in_primer", 'pool', 'time', 'age_group', 'age_treatment_time_mock_together'
  ), 
  colVarNums = c(), 
  save = 'ALI.CTCA'
)

print(plot_ranef(
  result$ranef,
  vars = list(age_treatment_time_mock_together = c('Paediatric-mock', 'Paediatric-SARS-4h', 'Paediatric-SARS-24h', 'Paediatric-SARS-72h')), 
  references = list(age_treatment_time_mock_together = 'Paediatric-mock'),
  celltype_order = rev_order, 
  maxFC=3, 
) + theme(axis.text.x=element_text(angle=45,hjust=1))
+ theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + ggtitle("Paediatric"))

print(plot_sdse(result$sdse, 'souporcell_cluster'))


# adult --------------------------------------------------------------
result <- CellTypeCompositionAnalysis(
  meta_for_GLMM_adult, 
  'souporcell_cluster', 
  'annotation_v5', 
  colVarCats=c(
    'donor_id', 'treatment', "kit_version", "spike.in_primer", 'pool', 'time', 'age_treatment_time_mock_together'
  ), 
  colVarNums = c(), 
  save = 'ALI.CTCA'
)

print(plot_ranef(
  result$ranef,
  vars = list(age_treatment_time_mock_together = c('Adult-mock', 'Adult-SARS-4h', 'Adult-SARS-24h', 'Adult-SARS-72h')), 
  references = list(age_treatment_time_mock_together = 'Adult-mock'),
  celltype_order = rev_order, 
  maxFC=3, 
) + theme(axis.text.x=element_text(angle=45,hjust=1))
+ theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + ggtitle("Adult"))

print(plot_sdse(result$sdse, 'souporcell_cluster'))

# elderly --------------------------------------------------------------
result <- CellTypeCompositionAnalysis(
  meta_for_GLMM_elderly, 
  'souporcell_cluster', 
  'annotation_v5', 
  colVarCats=c(
    'donor_id', 'pool', 'age_treatment_time_mock_together'
  ), 
  colVarNums = c(), 
  save = 'ALI.CTCA'
)

print(plot_ranef(
  result$ranef,
  vars = list(age_treatment_time_mock_together = c('Elderly-mock', 'Elderly-SARS-4h', 'Elderly-SARS-24h', 'Elderly-SARS-72h')), 
  references = list(age_treatment_time_mock_together = 'Elderly-mock'),
  celltype_order = rev_order, 
  maxFC=3, 
) + theme(axis.text.x=element_text(angle=45,hjust=1))
+ theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + ggtitle("Elderly"))

print(plot_sdse(result$sdse, 'souporcell_cluster'))


# Example from Ni ---------------------------------------------------------
result <- CellTypeCompositionAnalysis(
  obs_tbl1,
  'sample_id',
  'v6_annot2',
  colVarCats=c(
    'donor', 'Sex', 'Age_bin', 'Smoker', 'Kit_version', 'Sample_location', 'dataset'
  ),
  colVarNums=c(),
  extra_term='COVID_status:Group:Celltype',
  save='pooled_covid_airway.CTCA'
)

print(plot_ranef(
  result$ranef,
  vars=list(
    Sample_location=c('Nose', 'Trachea', 'Bronchi'),
    Age_bin=c('Neonate', 'Infant', 'Young child', 'Child', 'Adolescent', 'Adult', 'Eldly'),
    `COVID_status:Group`=c(
      'Healthy,Adult',
      'COVID+,Adult',
      'Post-COVID,Adult',
      'Healthy,Ped',
      'COVID+,Ped',
      'Post-COVID,Ped'
    )
  ),
  maxFC=3
) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))