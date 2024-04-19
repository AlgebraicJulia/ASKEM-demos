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
