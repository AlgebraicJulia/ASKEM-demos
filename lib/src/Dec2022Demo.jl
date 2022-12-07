module Dec2022Demo

export formSIRD, formTVParams, solveODE, zeroVal, runControlOptim, makeK, runControlAuto, draw, sig, invsig,
    formInfType, augLabelledPetriNet, sirdAugStates, typeSIRD, makeMultiAge, typeAge, 
    typed_stratify, formVax, vaxAugStates, typeVax, writeMdlStrat

# Common Imports
using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms

using AlgebraicPetri

# Imports for S2
using Catlab.Present
using AlgebraicPetri: Graph

# Imports for S1
using ModelingToolkit
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks
        
# using Random
using OrdinaryDiffEq
using Optimization
using OptimizationOptimisers
using LinearAlgebra    

import ..Oct2022Demo: sumvarsbyname, MakeReactionSystem
# using .Upstream: vectorfield
    
#*******************
# Functions for S1 *
#*******************

function formSIRD()
    SIRD = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	)
    return SIRD
end

sig(x) = 1/(1+exp(-x))
invsig(x) = log(x/(1-x))

function formTVParams(p,alpha,t_start)
    new_p_i(u,t) = p[1] * (1 - sig(alpha)*sig(t-t_start)) 
    new_p_r(u,t) = p[2]
    new_p_d(u,t) = p[3]
    p_t = [new_p_i,new_p_r,new_p_d]
    return p_t
end

function solveODE(tv_prob)
    tv_sol = solve(tv_prob, Tsit5())
end

function zeroVal()
    return 0
end

function runControlOptim(SIRD, tv_prob, p, tspan, hosp_rt, thresh_H, t_start, alpha_init)
    t = tspan[1]:tspan[2]
    
    function obsHfromI(model::AbstractLabelledPetriNet, sol, sample_times, hosp_rt)
        hosp_sample_vals = hosp_rt * sumvarsbyname(model, :I, sol, sample_times)
        labels = ["H"]
        return hosp_sample_vals, labels
    end
    
    function predictH(alpha)
        new_p_t = formTVParams(p,alpha,t_start)
        tv_sol = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5(), tstops=t)    
        vals, _ = obsHfromI(SIRD, tv_sol, t[1:findlast(t .<= t[end])], hosp_rt)
        return vals'
    end

    function l_func(control_params)
        pred_H = predictH(control_params[1])
        breach = sum(pred_H .> thresh_H)/length(pred_H)
        l = sig(control_params[1]) + exp(20*breach) - 1
        return l
    end

    optf = Optimization.OptimizationFunction((x, p) -> l_func(x), Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optf, alpha_init) 
    result1 = Optimization.solve(optprob, ADAM(0.1),  maxiters = 300) 
    # optprob2 = remake(optprob,u0 = result1.u)
    # result2 = Optimization.solve(optprob2, Optim.BFGS(initial_stepnorm=0.01)) 

    alpha = result1.u[1]
    new_p_t = formTVParams(p,alpha,t_start)
    tv_sol_cntrl = solve(remake(tv_prob,tspan=tspan,p=new_p_t), Tsit5())
    obs_hosp = tv_sol_cntrl(t)[2,:]*hosp_rt

    return alpha, tv_sol_cntrl, t, obs_hosp
end


function makeK(u,p)
    A = [1-p[1]*u[2] -p[1]*u[1] 0 0 2*u[1]u[2]*p[1] 0 0 0;
      p[1]*u[2] 1+p[1]*u[1]-p[2]-p[3] 0 0 0 -2*u[1]*u[2]*p[1]+u[2]*p[2]+u[2]*p[3] 0 0;
      0 p[2] 1 0 0 0 -p[2]*u[2] 0;
      0 p[3] 0 1 0 0 0 -p[3]*u[2];
      0 0 0 0 1 0 0 0;
      0 0 0 0 0 1 0 0;
      0 0 0 0 0 0 1 0;
      0 0 0 0 0 0 0 1]
  
    B = [-u[2]*u[3] 0 0;
    u[2]*u[3] -u[2] -u[2];
    0 u[2] 0;
    0 0 u[2];
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0]
  
    Q = zeros(8,8)
    Q[1,1] = 1
    Q[2,2] = 1
    Q[4,4] = 1
  
    R = Matrix{Float64}(I,3,3)
   
    Qf = Q
  
    PN = Q + A'*Qf*A - A'*Qf*B*inv(R+B'*Qf*B)*B'*Qf*A
    return K = inv(R + B'*PN*B)*B'*PN*A  
  end

function runControlAuto(rxn,u0,p,tspan)
    u_prev = u0
    p_prev = p
    t = tspan[1]:tspan[2]
    t_prev = t[1]
    sol_discr = zeros(length(u_prev),length(t))
    sol_discr[:,1] = u_prev
    for ii in 2:length(t)
        t_curr = t[ii]
        K_curr = makeK(u_prev,p_prev)
        u_aug = [u_prev; [1,1,1,1]]
        p_curr = -K_curr*u_aug
        prob = ODEProblem(rxn, u_prev, (t_prev,t_curr), p_curr)
        sol_curr = solve(prob, Tsit5())
        u_curr = sol_curr(t_curr)
        sol_discr[:,ii] = u_curr
        u_prev = u_curr
        p_prev = p_curr
        t_prev = t_curr
    end
    return sol_discr, t

end

#**********************
# Draw wiring diagram *
#**********************

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

#*******************
# Functions for S2 *
#*******************

# Type system model
types′ = LabelledPetriNet([:Pop],
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop))
types = map(types′, Name=name->nothing)

# Load unnamed type system
function formInfType()
    return types
end
  
# Parts of type system for ease of reference
s, = parts(types′, :S)
t_interact, t_disease, t_strata = parts(types′, :T)
i_interact1, i_interact2, i_disease, i_strata = parts(types′, :I)
o_interact1, o_interact2, o_disease, o_strata = parts(types′, :O);

# SIRD model already defined in S1
#=
function formSIRD()
    SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	)
    return SIRD_aug
end
=#

# Augment a LabelledPetriNet model with identity transitions for a set of states
function augLabelledPetriNet(lpn,states_to_aug)
  state_idx = Dict()
  for ii in 1:ns(lpn)
      state_idx[sname(lpn,ii)] = ii
  end
  for curr_s in lpn[:sname]
    if curr_s in states_to_aug
      i = add_transition!(lpn, tname=:id)
      add_inputs!(lpn, 1, i, state_idx[curr_s])
      add_outputs!(lpn, 1, i, state_idx[curr_s])
    end
  end
  return lpn
end

# Set of states to augment in SIRD
function sirdAugStates()
    return [:S, :I, :R]
end


# Typed SIRD model
function typeSIRD(SIRD_aug, types)
    SIRD_aug_typed = ACSetTransformation(SIRD_aug, types,
    S = [s, s, s, s],
    T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata],
    I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata],
    O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata],
    Name = name -> nothing 
    )
    @assert is_natural(SIRD_aug_typed)
    return SIRD_aug_typed
end

# Assemble a Multi-Age augmented model
function makeMultiAge(n;f_aug=true)
    lstates = []
    ltrans = []
    for ii in 1:n
        push!(lstates,Symbol("Age"*string(ii)))
    end
    for ii in 1:n
        for jj in 1:n
            push!(ltrans,Symbol("inf"*string(ii)*string(jj)) => ((lstates[ii],lstates[jj])=>(lstates[ii],lstates[jj])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("dis"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    if f_aug
        for ii in 1:n
            push!(ltrans,Symbol("id"*string(ii)) => ((lstates[ii])=>(lstates[ii])))
        end
    end
    MultiAge = LabelledPetriNet(lstates,ltrans...)
end

# Typed Multi-Age model
function typeAge(MultiAge_aug,types)
    num_ages = ns(MultiAge_aug)
    ma_S = repeat([s],num_ages)
    ma_T = repeat([t_interact],num_ages*num_ages)
    ma_I = repeat([i_interact1, i_interact2],num_ages*num_ages)
    ma_O = repeat([o_interact1, o_interact2],num_ages*num_ages)
    for ii in 1:num_ages
        push!(ma_T,t_disease)
        push!(ma_I,i_disease)
        push!(ma_O,o_disease)
    end
    for ii in 1:num_ages
        push!(ma_T,t_strata)
        push!(ma_I,i_strata)
        push!(ma_O,o_strata)
    end
    MA_aug_typed = ACSetTransformation(MultiAge_aug, types,
	  S = ma_S,
	  T = ma_T,
	  I = ma_I,
	  O = ma_O,
	  Name = name -> nothing 
    )
    @assert is_natural(MA_aug_typed)
    return MA_aug_typed
end

# Function to form stratified model
typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)


# Vaccination augmented model
function formVax()
    Vax_aug = LabelledPetriNet([:U, :V],
	  :infuu => ((:U, :U)=>(:U, :U)),
	  :infvu => ((:V, :U)=>(:V, :U)),
	  :infuv => ((:U, :V)=>(:U, :V)),
	  :infvv => ((:V, :V)=>(:V, :V)),
	  :vax => (:U => :V),
	)
    return Vax_aug
end 

# Set of states to augment in Vax
function vaxAugStates()
    return [:U, :V]
end

# Typed Vax model
function typeVax(Vax_aug,types)
    Vax_aug_typed = ACSetTransformation(Vax_aug, types,
    S = [s, s],
    T = [t_interact, t_interact, t_interact, t_interact, t_strata, t_disease, t_disease],
    I = [i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_interact1, i_interact2, i_strata, i_disease, i_disease],
    O = [o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_interact1, o_interact2, o_strata, o_disease, o_disease],
    Name = name -> nothing 
    )
    @assert is_natural(Vax_aug_typed)
    return Vax_aug_typed
end

# Write stratified model to file as JSON
function writeMdlStrat(mdl, file)
  write_json_acset(dom(mdl), file)
end



end