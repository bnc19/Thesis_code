
################################################################################
# This code extracts the 1000 non-vac simulations with lowest DENV2 incidence  #
# and highest DENV3 incidence -- it will only work if the full model is run    # 
# not the demo                                                                 #
################################################################################

setwd("Chapter4")

library(matrixStats)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
source("/R/process_simulations")

pathin = "outputs/" 
pathout = paste0(getwd(), "/outputs/processed/")

countries = c("BRA")
vacc.age = c(6)
vacc.cov = c(0.8) 
R0 = c("R01","R02","R03","R04","R05","R06","R07","R08","R09")
serop9 = c("sp9 10%", "sp9 20%", "sp9 30%", "sp9 40%", "sp9 50%", "sp9 60%", "sp9 70%", "sp9 80%", "sp9 90%")

# serotype specific disease in whole population 
v = c ("out_dis_sero_pop")

# number of simulations
nsim = 1000

# serotype to focus on 
serotype_name = 2

# "min" or "max" incidence 
inc = 'min'

folder_type = c("/VI_D15/" , "/VS_D15/")

# multipliers for dis and sdis
dis.scale = 1.2956
sdis.scale = c(0.6833, 1.0553)

new_folder = paste0(dis.scale, "_", sdis.scale[1], "_", sdis.scale[2], "/")

# iterate over folders 
out.list = list()

# path in and out 

for(h in 1: length(countries)){
  for(r in 1:length(R0)){

pathin.v = paste0(pathin, countries[h], "/vacc/")
pathin.nv = paste0(pathin, countries[h], "/no_vacc/")

pathin.vik = paste0(pathin.v, "vc_", vacc.cov, folder_type[1])
pathin.vsk = paste0(pathin.v, "vc_", vacc.cov, folder_type[2])
pathin.nvk = paste0(pathin.nv, "vc_", vacc.cov)

scenario = paste0(vacc.age,"_", R0[r]) # transmission settings 

# read in results 
result.nv = readRDS(paste0(pathin.nvk,"/no_vacc_res_", scenario,".rds"))
result.vi = readRDS(paste0(pathin.vik,"/vacc_res_", scenario,".rds"))
result.vs = readRDS(paste0(pathin.vsk,"/vacc_res_", scenario,".rds"))

# yearly serotype outcome (first 20 years)
t.vi = dis.scale*result.vi[[v]]
t.vs = dis.scale*result.vs[[v]]
t.nv = dis.scale*result.nv[[v]]

# get serotype cases across all scenarios 
all_cases = list(
nv_cases_1 = (t.nv[1:20,,]),
vs_cases_1 = (t.vs[1:20,,]),
vi_cases_1 = (t.vi[1:20,,]),

nv_cases_2 = (t.nv[21:40,,]),
vs_cases_2 = (t.vs[21:40,,]),
vi_cases_2 = (t.vi[21:40,,]),

nv_cases_3 = (t.nv[41:60,,]),
vs_cases_3 = (t.vs[41:60,,]),
vi_cases_3 = (t.vi[41:60,,]),

nv_cases_4 = (t.nv[61:80,,]),
vs_cases_4 = (t.vs[61:80,,]),
vi_cases_4 = (t.vi[61:80,,])
)


  # get cumulative serotype cases for each nv simulation 
  col_name = paste0("nv_cases_", serotype_name)   
  ind_col = which(names(all_cases) == col_name)
  cum = apply(all_cases[[ind_col]], c(2,3), sum)
  
  # get either the minimum or maximum incidence 
  if(inc  == "max"){
  # get the position of the nsim smallest cumulative incidences
  index = which(cum >= tail(sort(cum), nsim)[1], arr.ind = T)
  } else {
  # get the position of the nsim largest cumulative incidences 
    index = which(cum <= head(sort(cum), nsim)[nsim], arr.ind = T)
  }
  
  
  # create data frame labeling each indexed simulation using row number 
  ni_df = index %>% 
    as.data.frame() %>%
    rownames_to_column() 
  
  # create data frame of years (1:20) to match number of simulations 
  time_df = 
    data.frame(
      time = rep(1:20, nrow(ni_df)),
      rowname = as.character(rep(1:nrow(ni_df), each = 20))
    )
  
  # left join simulation indexing to years, so we can extract every year for simulations with highest D3 cases 
  ni_df_2 =  time_df %>% 
    left_join(ni_df) %>% 
    select(-rowname) %>% 
    as.matrix()
  
  # selected same indexing for all serotypes and simulations 
  ind_cases = lapply(all_cases, function(x) {
    x[ni_df_2]
  })
  


name = rep(c("nv", "ns", "ni"), 4)
serotypes = rep(paste0("DENV", 1:4), each = 3)

# add year, iteration number, simulation name and serotype
form_ind_cases = list()
for(i in 1:length(ind_cases)){
  form_ind_cases[[i]] = format_sero_dyn(ind_cases[[i]], name = name[i], serotypes=serotypes[i])
}

# single data frame out 
out = form_ind_cases %>%  
  bind_rows()

out.list[[paste0(countries[h], "_", R0[r])]] = out

}
}

saveRDS(out.list, paste0(pathout, "/serotype_dynamics_nsim_",nsim, "_", serotype_name, "_", inc, ".rds"))
