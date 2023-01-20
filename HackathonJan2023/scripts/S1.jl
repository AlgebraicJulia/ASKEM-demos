

#*****
# Q2 *
#*****

# SEIRD
function formSIRD()
    SIRD = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	)
    return SIRD
end

#*******
# Q3 *
#*******

# CHIMESVIIvR

# SEIRD-V

#***
# Q3b
#***
# Pairwise comparison

# Three-way comparison

# Plot structural comparisons

#***
# Q3c
#***
# Sim all three models
# i)  Sim w/ vax effic 75%, pop vax 10%
# ii) Sim w/ vax effic 75%, pop vax 10%
# Compare sim outputs

#*****
# Q4 *
#*****
# Create eq-wt form_ensemble_model
# Sim for 3ci and 3cii
# Compare w/ indiv run_model_selection

#*****
# Q5 *
#*****
# Sensitivity sensitivity analysis

#*****
# Q6 *
#*****
# Stratify one model (all models) by age.

#***
# Q6d
#***
# Sim model(s) for...
# i)   High vax rt (+80%) for 65+yo and low vax rt (<15%) for all others
# ii)  High vax rt for all grps
# iii) Repeat i and ii with 20% decrease in contact for school-age children
# iv)  Compare outputs

