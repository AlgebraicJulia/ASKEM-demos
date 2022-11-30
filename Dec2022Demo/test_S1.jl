types′ = LabelledPetriNet([:Pop],
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop))
types = map(types′, Name=name->nothing)


SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));
AlgebraicPetri.Graph(SIRD)

SIRD_typed = homomorphism(SIRD, types;
initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
type_components=(Name=x->nothing,))
@assert is_natural(SIRD_typed)


p_i(u,t) = 0.000025-0.000025*cos(t/100*2*pi)
# sigmoid(x) = 1/(1-exp(-x))
# p_i(u,t) = 0.00005*(1/(1+exp(-(t-30))))
p_r(u,t) = 0.005+0.005*cos(t/100*4*pi+pi/4)
p_d(u,t) = 0.001+0.001*cos(t/100*8*pi+pi/2)

p_t = [p_i,p_r,p_d]

tv_rxn = vectorfield(SIRD)
tspan = (0,1000)
u0 = [999,1,0,0]
tv_prob = ODEProblem(tv_rxn, u0, tspan, p_t)
tv_sol = solve(tv_prob, Tsit5())
plot(tv_sol)

p_i(u,t) = 0.000025-0.000025*cos(t/100*2*pi)
p_r(u,t) = 0.01
p_d(u,t) = 0.002
p_t = [p_i,p_r,p_d]


    
#****

function evalparam(β::Function,t)
  β(t)
end
function evalparam(β::Number,t)
  β
end

rates[i] => evalparam(rates[i],t)


rxn = MakeReactionSystem(model)

rxn = vectorfield(model) # f!(du, u, p, t)

tspan = (sample_times[1], sample_times[end])

prob = ODEProblem(rxn, u0, tspan, p_init)

p_estimate = p_init
loss = 0
sol_estimate = Nothing

AtoB = LabelledPetriNet([:A, :B], :recover => (:A=>:B));
p_r(t) = 0.05+0.01*cos(t/100*4+pi/4)

function p(t) 
  [p_r(t)]
end

p_t = Dict()
p_t[:inf] = p_i
p_t[:recover] = p_r
p_t[:death] = p_d
