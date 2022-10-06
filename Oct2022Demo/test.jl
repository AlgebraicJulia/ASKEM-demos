include(Oct2022Demo.jl)

infectious_type = LabelledPetriNet([:Pop],
:interact=>((:Pop, :Pop)=>(:Pop, :Pop)), 
:t_disease=>(:Pop=>:Pop),
:t_strata=>(:Pop=>:Pop)
)
s, = parts(infectious_type, :S)
t_interact, t_disease, t_strata = parts(infectious_type, :T)
i_interact1, i_interact2, i_disease, i_strata = parts(infectious_type, :I)
o_interact1, o_interact2, o_disease, o_strata = parts(infectious_type, :O);
infectious_type = map(infectious_type, Name=name->nothing); 

SIR = LabelledPetriNet([:S, :I, :R],
:inf => ((:S, :I)=>(:I, :I)),
:rec => (:I=>:R),
:id => (:S => :S),
:id => (:I => :I),
:id => (:R => :R)
)
typed_SIR = ACSetTransformation(SIR, infectious_type,
S = [s, s, s],
T = [t_interact, t_disease, t_strata, t_strata, t_strata],
I = [i_interact1, i_interact2, i_disease, i_strata, i_strata, i_strata],
O = [o_interact1, o_interact2, o_disease, o_strata, o_strata, o_strata],
Name = name -> nothing # specify the mapping for the loose ACSet transform
);
@assert is_natural(typed_SIR)

quarantine = LabelledPetriNet([:Q, :NotQ],
:id => (:Q => :Q),
:id => (:NotQ => :NotQ),
:enter_q => (:NotQ => :Q),
:exit_q => (:Q => :NotQ),
:interact => ((:NotQ, :NotQ) => (:NotQ, :NotQ))
)
typed_quarantine = ACSetTransformation(quarantine, infectious_type,
S = [s, s],
T = [t_disease, t_disease, t_strata, t_strata, t_interact],
I = [i_disease, i_disease, i_strata, i_strata, i_interact1, i_interact2],
O = [o_disease, o_disease, o_strata, o_strata, o_interact1, o_interact2],
Name = name -> nothing
)
@assert is_natural(typed_quarantine)

true_mdl = typed_stratify(typed_SIR,typed_quarantine)
true_p = vcat(repeat([1e-4], 14), repeat([0.01], 6))
true_p = vcat([1.1e-4,1.3e-4,1.2e-4,1.4e-4,1.5e-4,1.2e-4,
             1.e-4,1.3e-4,1.5e-4,1.2e-4,1.7e-4,1.6e-4,
             1.2e-4,1.1e-4,], [0.01,0.02,0.01,0.01,0.01,0.02])
u0 = [0.0,0.0,0.0,0.0,0.0,999.0,1.0,0.0,0.0,0.0] 
tspan = (0.0,250.0)

true_sol = simulate(true_mdl,u0,true_p,tspan)

sample_data, sample_times, prob_real, sol_real = generate_data(true_mdl, true_p, u0, tspan, 50, obsSIRD)

plt = plot(sol_real, lw=2, label=reshape(map(string, true_mdl[:, :sname]), 1, ns(true_mdl)))
plot!(sample_times, sample_data, seriestype=:scatter, label="")
