module Mira 
export load_mira, load_mira_curated

function load_mira(fname)
  mdl_orig_mira = read_json_acset(OrigMIRANet{Any,Any,Any,Any,Any,Any}, fname)
  mdl_orig_mira[:sname] .= Symbol.(mdl_orig_mira[:sname])
  mdl_orig_mira[:tname] .= Symbol.(mdl_orig_mira[:tname])
  mdl_orig = LabelledPetriNet(mdl_orig_mira);
  return mdl_orig
end

function load_mira_curated(fname, fillrate=1.0)
  mdl_disease_mira = read_json_acset(MIRANet{Any,Any,Any,Any,Any,Any}, fname)
  input_rates = deepcopy(mdl_disease_mira[:rate])
	min_rate = minimum(input_rates[map(!isnothing,input_rates)])
  filled_rates = input_rates
  filled_rates[map(isnothing,input_rates)] = repeat([fillrate],sum(map(isnothing,input_rates)))
  mdl_disease_mira[:rate] = filled_rates
  mdl_disease_mira[:concentration] = 0.0
  mdl_disease_mira[:sname] .= Symbol.(mdl_disease_mira[:sname])
  mdl_disease_mira[:tname] .= Symbol.(mdl_disease_mira[:tname])
  return mdl_disease_mira
end


end # module 
