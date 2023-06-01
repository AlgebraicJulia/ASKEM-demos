using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms


draw(d::WiringDiagram) = to_graphviz(d,
    orientation=LeftToRight,
    labels=true, label_attr=:xlabel,
    node_attrs=Graphviz.Attributes(
      :fontname => "Courier",
    ),
    edge_attrs=Graphviz.Attributes(
      :fontname => "Courier",
    )
)

# Form Workflow presentation of FreeBiproductCategory
@present Workflow(FreeBiproductCategory) begin
    (File,LRN,LPN,TypedLPN,StrataSpec,LPNss,ObsFunc,ParamVec,StateVec,TSpan,NumTS,SampleData,SampleTimes,ODEProb,ODESol,Labels,Loss)::Ob 
    LoadLRN::Hom(File,LRN)
    
    Homomorph::Hom(LPN⊗LPN,TypedLPN)
    StrataSpecify::Hom(TypedLPN⊗StrataSpec,LPNss)
    Stratify::Hom(LPNss⊗LPNss⊗LPN,LPN⊗ObsFunc) 

    GenerateData::Hom(LPN⊗ObsFunc⊗ParamVec⊗StateVec⊗TSpan⊗NumTS,SampleData⊗SampleTimes⊗ODEProb⊗ODESol⊗SampleData⊗Labels)
    Calibrate::Hom(LPN⊗ObsFunc⊗StateVec⊗ParamVec⊗SampleData⊗SampleTimes,ParamVec⊗ODESol⊗Loss)

    # FormControl
end

# Form wiring diagram of load_stratify_calibrate_control Workflow
load_stratify_calibrate_control = @program Workflow (f_disease::File,f_strat::File,f_type::File,
                                                        spec_disease::StrataSpec,spec_strat::StrataSpec,
                                                        true_p::ParamVec,true_u0::StateVec,tspan::TSpan,num_ts::NumTS,
                                                        p_init::ParamVec) begin # 
    mdl_disease = LoadLRN(f_disease)
    mdl_strat = LoadLRN(f_strat)
    mdl_type = LoadLRN(f_type)

    # Form stratified model
    mdl_disease_typed = Homomorph(mdl_disease,mdl_type)
    mdl_strat_typed = Homomorph(mdl_strat,mdl_type)
    disease_ss = StrataSpecify(mdl_disease_typed,spec_disease)
    strat_ss = StrataSpecify(mdl_strat_typed,spec_strat)
    mdl_stratified, obs_func = Stratify(disease_ss,strat_ss,mdl_type) 

    # Simulate data
    sample_data, sample_times, prob_true, sol_true, noiseless_data, data_labels = GenerateData(mdl_stratified, obs_func, true_p, true_u0, tspan, num_ts)

    # Calibrate
    p_est, sol_est, loss = Calibrate(mdl_stratified, obs_func, true_u0, p_init, sample_data, sample_times)

    # Form controler

    return  p_est, sol_est, loss, mdl_stratified, obs_func, sample_data, sample_times, prob_true, sol_true, noiseless_data, data_labels
end

# Display wiring diagram of workflow
draw(load_stratify_calibrate_control)

# Write diagram to file as JSON
write_json_acset(load_stratify_calibrate_control.diagram, "diagram_load_strat_calib_cntrl.json")


