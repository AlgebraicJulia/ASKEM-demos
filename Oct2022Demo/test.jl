include("Oct2022Demo.jl")

types′ = LabelledPetriNet([:Pop],
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop))
types = map(types′, Name=name->nothing)

SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D))

SIRD_typed = homomorphism(SIRD, types;
    initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
    type_components=(Name=x->nothing,))

@assert is_natural(SIRD_typed)

SIRD_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])

Quarantine = LabelledPetriNet([:Q,:NQ],
    :quarantine => ((:NQ)=>(:Q)),
    :unquarantine => ((:Q)=>(:NQ)))

Quarantine_typed = homomorphism(Quarantine, types;
    initial=(T=[3,3],), type_components=(Name=x->nothing,))

@assert is_natural(Quarantine_typed)

Quarantine_ss = StrataSpec(Quarantine_typed, [[:disease], [:disease,:infect]])

true_mdl = stratify(SIRD_ss, Quarantine_ss, types′)
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
