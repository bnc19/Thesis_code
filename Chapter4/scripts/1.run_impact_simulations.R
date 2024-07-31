
########### Set up to run on the high performance computing cluster ############
# library(hipercow)
# hipercow_configure(driver = "windows")
# packages = c("dust", "odin", "odin.dust", "denguetak")
# my_sources = c("R/run_who_cluster.r","R/process_demog.r","R/config_params.r")
# hipercow_provision()
# hipercow_environment_create(sources = my_sources,
#                             packages = packages)
################################################################################ 


############################ Set up to run locally #############################
# # install packages
# install.packages("devtools")
# install.packages("pkgbuild") 
# install.packages("pkgload")
# 
#devtools::install_github("mrc-ide/dust")
#devtools::install_github("mrc-ide/odin")
#devtools::install_github("mrc-ide/odin.dust")
#devtools::install_github("NeilFerguson/denguetak")
################################################################################ 

setwd("Chapter4")

library(odin.dust)
library(dust)
library(odin)

# source files  
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

################################################################################ 
################################### DEMO CODE ##################################
################################################################################ 

# To reconstruct the results in the thesis, I drew 200 samples from the 
# posterior distribution of the parameters of the VE model and ran 50 stochastic 
# simulations for each such posterior sample, giving 10,000 simulations per
# scenario in total. 

# The scenarios run include:

#  - Vaccine mechanism of protection: only disease or also against infection (12, 3, 12_3)
#  - Duration of vaccine decay: 15 years or 5 years
#  - Age of vaccination: 6 to 12 years
#  - Transmission intensity setting (average 9 years seroprevalence): 10% to 90%
#  - Vaccine coverage: 20%, 40%, 60%, 80%
#  - Pre-vaccination screening to identify seropositive individuals: yes/no 

# The code in this repo shows how to run four scenarios (Vaccine 
# protects disease only for 15 years, vaccination of age 6 with 80% coverage,
# without pre-vaccination screening, in transmission settings with 20%, 40%, 60%,
# and 80% 9-year-old seropositivity).

# For these scenarios, the demo samples the posterior VE estimates 50 times, 
# and run 10 stochastic simulations per sample.


# The code to run the full model across all scenarios on a HPC is below the demo.
# Each scenario is run as a separate job, taking ~12 hours to run 200 posterior 
# samples and 50 simulations per sample. 

# model set up
ns = 50 # number of posterior samples of VE estimates to use
np = 10 # number of stochastic simulations per posterior sample
nt = 50 # number of threads to run at once

pf = "post.samp.csv"  # VE param estimates posterior file
outf = "demo"     # location of results 

# create output folders 
dir.create("demo", showWarnings = F)
dir.create("demo/BRA", showWarnings = F)
dir.create("demo/BRA/vacc", showWarnings = F)
dir.create("demo/BRA/no_vacc", showWarnings = F)
dir.create("demo/BRA/equil", showWarnings = F)

##### scenario parameters 

i_country = 1 # country index (1 = Brazil, 2 = Phillipines)
i_R0 = c(2,4,6,8)  # transmission strength index (1:9, 10:90% average 9-year-old seropositivty)
vacc.age = 6 # vaccination age 
vacc_cov = 0.8 # vaccination coverae 
vacc_by_serostatus = 0 # pre-vaccination screening (0 = no, 1 = yes)
DoVEInf = 0 # vaccine efficacy against only disease (0) or also against infection (1)
mdy = 15 # duration of vaccine efficacy decay (years)

# 1) Equilibriate the model for 150 years assuming static demography prior 1950 
#    and non-stationary demography thereafter (~ 1 hour 50 minute run time)


for(i in i_R0){
run_WHO_eq(i_country = i_country,
           i_R0 = i,
           ns = ns,
           nparticles = np,
           posterior_file = pf,
           nthreads = nt,
           outf = outf)
}


# 2) Run the counterfactual, no vaccination scenario (run time ~ 12 minutes)

for(i in i_R0){
run_WHO_novacc(
  i_country = i_country,
  i_R0 = i,
  ns = ns,
  nparticles = np,
  posterior_file = pf,
  nthreads = nt,
  vacc.age = vacc.age,
  vacc_cov = vacc_cov,
  vacc_by_serostatus = vacc_by_serostatus,
  outf = outf
)}



# 3) Run the vaccine scenario (run time ~12 minutes)

for(i in i_R0){ 
run_WHO_vacc(
  i_country = i_country,
  i_R0 = i,
  ns = ns,
  nparticles = np,
  posterior_file = pf,
  nthreads = nt,
  vacc.age = vacc.age,
  vacc_cov = vacc_cov,
  DoVEInf = DoVEInf,
  vacc_by_serostatus = vacc_by_serostatus,
  mdy = mdy,
  outf = outf
)
}


# 
# ################################################################################ 
# ######### CODE TO RECREATE RESULTS PRESENTED IN THE THESIS  ####################
# ################################################################################ 
# 
# dir.create("outputs", showWarnings = F)
# dir.create("outputs/BRA", showWarnings = F)
# dir.create("outputs/BRA/vacc", showWarnings = F)
# dir.create("outputs/BRA/no_vacc", showWarnings = F)
# dir.create("outputs/BRA/equil", showWarnings = F)
# dir.create("outputs_vbs", showWarnings = F)
# dir.create("outputs_vbs/BRA", showWarnings = F)
# dir.create("outputs_vbs/BRA/vacc", showWarnings = F)
# dir.create("outputs_vbs/BRA/no_vacc", showWarnings = F)
# dir.create("outputs_vbs/BRA/equil", showWarnings = F)
# dir.create("outputs_3VI", showWarnings = F)
# dir.create("outputs_3VI/BRA", showWarnings = F)
# dir.create("outputs_3VI/BRA/vacc", showWarnings = F)
# dir.create("outputs_3VI/BRA/no_vacc", showWarnings = F)
# dir.create("outputs_3VI/BRA/equil", showWarnings = F)
# dir.create("outputs_12VI", showWarnings = F)
# dir.create("outputs_12VI/BRA", showWarnings = F)
# dir.create("outputs_12VI/BRA/vacc", showWarnings = F)
# dir.create("outputs_12VI/BRA/no_vacc", showWarnings = F)
# dir.create("outputs_12VI/BRA/equil", showWarnings = F)

# ns = 200
# np = 50
# nt = 50
# pf = "post.samp.csv"
# outf = "Z:/outputsBCD/impact_BCD/outputs"
# outf2 = "Z:/outputsBCD/impact_BCD/outputs"
# outf = outf2
# 
# ## vacc runs
# 
# job = list()
# k = 0
# vbs = 0 
# 
# # age for VS D15 0.8
# for(vc in c(0.8)) {   
#   for(va in c(6:12)) { # make 6 - 12 if want age differences 
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(0)) {
#         for(mdy in c(15)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
# 
# 
# # VS D5 
# 
# job = list()
# k = 0
# vbs = 0 
# 
# for(vc in c(0.8)) {   
#   for(va in c(6)) {
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(0)) {
#         for(mdy in c(5)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
# 
# 
# #  VI coverage for age 6, decay 15
# job = list()
# k = 0
# vbs = 0 
# 
# for(vc in c(0.2,0.4, 0.6, 0.8)) {   
#   for(va in c(6)) { 
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(1)) {
#         for(mdy in c(15)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
# 
# # SA on VI 3
# 
# outf = "Z:/outputsBCD/impact_BCD/outputs_3VI"
# outf2 = "Z:/outputsBCD/impact_BCD/outputs_3VI"
# 
# 
# job = list()
# k = 0
# vbs = 0 
# 
# for(vc in c(0.8)) {   
#   for(va in c(6)) {  
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(1)) {
#         for(mdy in c(15)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
# 
# 
# 
# 
# # SA on VI 12
# 
# outf = "Z:/outputsBCD/impact_BCD/outputs_12VI"
# outf2 = "Z:/outputsBCD/impact_BCD/outputs_12VI"
# 
# 
# job = list()
# k = 0
# vbs = 0 
# 
# for(vc in c(0.8)) {   
#   for(va in c(6)) {  
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(1)) {
#         for(mdy in c(15)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
# 
# 
# ## RUN SCREENING 
# outf = "Z:/outputsBCD/impact_BCD/outputs_vbs"
# outf2 = "Z:/outputsBCD/impact_BCD/outputs_vbs"
# 
# 
# 
# job = list()
# k = 0
# vbs = 1
# 
# # age for VS D15 0.8
# for(vc in c(0.8)) {   
#   for(va in c(6)) { # make 6 - 12 if want age differences 
#     ## no vacc runs
#     for(i in 1:9) {
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(1,i,ns,np,pf,nt,va,vc,vbs,outf))
#       k = k+1
#       job[[k]] = task_create_expr(run_WHO_novacc(2,i,ns,np,pf,nt,va,vc,vbs,outf))
#       ## vacc runs
#       for(vei in c(0)) {
#         for(mdy in c(15)) {
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(1,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#           k = k+1
#           job[[k]] = task_create_expr(run_WHO_vacc(2,i,ns,np,pf,nt,va,vc,vei,vbs,mdy,outf))
#         }
#       }
#     }
#   }
# }
