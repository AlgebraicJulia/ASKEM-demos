using Base.Threads
using Catlab, Catlab.Theories
using Catlab.CategoricalAlgebra
using Catlab.Graphics
using Catlab.Graphics: Graphviz
import Catlab.CategoricalAlgebra: migrate!
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Programs.RelationalPrograms
using Catlab.Syntax
using Catlab.Programs.ParseJuliaPrograms

using Catlab.Present
using AlgebraicPetri
using AlgebraicPetri: Graph
using AlgebraicPetri.Epidemiology
using AlgebraicPetri.BilayerNetworks

import ASKEM.Upstream: presentationToLabelledPetriNet
import Catlab.WiringDiagrams.DirectedWiringDiagrams: WiringDiagramACSet
using Catlab.CategoricalAlgebra.CSets

using ModelingToolkit
using OrdinaryDiffEq
using DifferentialEquations

using JSON, CSV, DataFrames

### construct a simple concurrent interpreter from a program/wiring diagram

function interpreter(prog, env)
    ws = wires(prog)
    bs = boxes(prog)
    channels = [ Channel(0) for _ in ws ]
    function put_value!(box, port, value)
        # put value into every wire from the given source
        for i in 1:length(ws)
            w = ws[i]
            if w.source.box == box && w.source.port == port
                put!(channels[i], value)
            end
        end
    end
    n_inputs = length(input_ports(prog))
    n_outputs = length(output_ports(prog))
    set_input!(port, value) = put_value!(-2, port, value)
    set_inputs!(values...) = foreach(i->set_input!(i, values[i]), 1:length(values::NTuple{n_inputs,Any}))
    function channels_for_box(b, kind)
        if b == -1
            n = n_outputs
        elseif b == -2
            n = n_inputs
        else
            n = length(kind === :target ? bs[b].input_ports : bs[b].output_ports)
        end
        chans = Vector{Channel{Any}}(undef, n)
        for j in 1:length(ws)
            w = ws[j]
            if getfield(w, kind).box == b
                chans[getfield(w, kind).port] = channels[j]
            end
        end
        @assert all(isassigned(chans, i) for i in 1:length(chans))
        return chans
    end
    for i in 1:length(bs)
        n_box_outputs = length(bs[i].output_ports)
        input_chans = channels_for_box(i, :target)
        workfun = bs[i].value
        function execute(inputs)
            #=
            str = "calling $workfun"
            if !isempty(inputs)
                str *= " on " * join(inputs, ", ")
            end
            sleep(rand()*0.1)
            println(str); flush(stdout)
            #for p in 1:n_box_outputs
            #    put_value!(i, p, "box $i result $p")
            #end
            =#
            try
                put_value!(i, 1, env[workfun](inputs...))
            catch e
                if e isa InvalidStateException || e isa InterruptException
                else
                    println(e)
                end
                rethrow()
            end
        end
        @spawn while true
            inputs = Any[ take!(c) for c in input_chans ]
            execute(inputs)
        end
    end
    out_chans = channels_for_box(-1, :target)
    get_outputs!() = map(take!, out_chans)
    stop!() = foreach(close, channels)
    (; set_inputs!, get_outputs!, stop!)
end

# example:
# sim = interpreter(prog)
# sim.set_inputs!(1,2,3,4)
# sim.get_outputs!()

### compiling to MTK

make_depvar(p,t) = :($p($t))

function compile_mtk(bn::Union{AbstractLabelledBilayerNetwork, AbstractBilayerNetwork})
  varstmt = :(@variables t)
  varnames = bn[:variable]
  append!(varstmt.args, make_depvar.(bn[:variable], :t))

  paramstmt = :(@parameters)
  params = bn[:parameter]
  append!(paramstmt.args, bn[:parameter])

  diffstmt = :(_D_ = Differential(t))

  ϕs = map(parts(bn, :Box)) do b
    vars = map(incident(bn, b,:call)) do i
      j = bn[i, :arg]
      return bn[j, :variable]
    end
    p = :(*($(bn[b, :parameter])))
    append!(p.args, vars)
    return :($(Symbol("ϕ$b")) = $p)
  end

  infs = map(parts(bn, :Qout)) do tv
    vars = map(incident(bn, tv, :infusion)) do wa
      j = bn[wa, :influx]
      return Symbol("ϕ$j")
    end
    p = :(+())
    append!(p.args, vars)

    # same for the outfluxes
    vars = map(incident(bn, tv, :effusion)) do wn
      j = bn[wn, :efflux]
      return :(- $(Symbol("ϕ$j")))
    end
    append!(p.args, vars)
    return p
  end

  zparts = zip(bn[:tanvar], infs)

  eqns = [:(_D_($tanvar) ~ $rhs) for (tanvar, rhs) in zparts]
  eq = :([])
  append!(eq.args, eqns)
  eqnstmt = :(eqs = $eq)

  varnameexpr = Expr(:tuple, varnames...)
  parnameexpr = Expr(:tuple, params...)

  return quote
    $varstmt
    $paramstmt
    $diffstmt
    $(ϕs...)
    $eqnstmt
    return $varnameexpr, $parnameexpr, ODESystem(eqs, t, name=:PetriNet)
  end
end

### diagram/model I/O

deserialize_wiringdiagram(filepath::String) = deserialize_wiringdiagram!(read_json_acset(WiringDiagramACSet{Symbol,Any,Any,Any}, filepath))
deserialize_wiringdiagram!(dwd) = begin
  convsymbol(dwd, key) = begin
    dwd[key] .= Symbol.(dwd[key])
  end
  
  dwd[:box_type] .= Box{Symbol}
  convsymbol(dwd, :in_port_type)
  convsymbol(dwd, :out_port_type)
  convsymbol(dwd, :outer_in_port_type)
  convsymbol(dwd, :outer_out_port_type)
  dwd[:value] .= map(Symbol,dwd[:value])
  wd_acset2 = WiringDiagramACSet{Symbol,Any,Any,DataType}()
  copy_parts!(wd_acset2,dwd)
  return WiringDiagram{ThBiproductCategory, Symbol, Any, Any}(wd_acset2, :read)
end

### simulation plan functional blocks

function formSIRD()
    SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
    return SIRD_aug
end

function json_to_mdl(json::String)
    parse_json_acset(LabelledPetriNet, json)
end

function mtk_simulate(model::LabelledPetriNet, states, params, timespan)
    bnsir = LabelledBilayerNetwork()
    migrate!(bnsir, model)
    mdl = eval(compile_mtk(bnsir))
    prob = ODEProblem(mdl[3], states, timespan, params)
    soln = solve(prob)
    return soln
end

function soln_to_json(soln::ODESolution)
    json(Dict(:time=>soln.t, :states=>soln.u))
end

function soln_to_csv(soln::ODESolution)
    io = IOBuffer()
    CSV.write(io, DataFrame(soln))
    String(take!(io))
end

@present ForecastWorkflow(FreeBiproductCategory) begin
    (ParamVec,StateVec,TSpan,String,Mdl,Soln)::Ob
    formSIRD::Hom(munit(),Mdl)
    mtk_simulate::Hom(Mdl⊗StateVec⊗ParamVec⊗TSpan,Soln)
    soln_to_csv::Hom(Soln,String)
end

#=
# baked in as forecast_plan.json
forecast = @program ForecastWorkflow (state::StateVec, params::ParamVec, ts::TSpan) begin
    mdl = formSIRD()
    soln = mtk_simulate(mdl, state, params, ts)
    return soln_to_csv(soln)
end

@program ForecastWorkflow (model::String, state::StateVec, params::ParamVec, ts::TSpan) begin
    mdl = json_to_mdl(model)
    soln = mtk_simulate(mdl, state, params, ts)
    return soln_to_csv(soln)
end
=#
