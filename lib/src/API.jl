using NetCDF

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

deserialize_wiringdiagram(json::String) = deserialize_wiringdiagram!(parse_json_acset(WiringDiagramACSet{Symbol,Any,Any,Any}, json))
function deserialize_wiringdiagram!(dwd)
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

function soln_to_netcdf(soln::ODESolution, filename)
    df = DataFrame(soln)
    time = df[!, "timestamp"]
    attrnames = filter(!=("timestamp"), names(df))
    data = Matrix(df[!, attrnames])
    nccreate(filename, "solution",
             "timesteps", time,
             "attributes", attrnames)
    ncwrite(data, filename, "solution")
    filename
end

function run_workflow(wiringdiagram, inputs, output_filename)
    wd = deserialize_wiringdiagram(wiringdiagram)
    args = JSON.parsefile(inputs)
    sim = interpreter(wd, (; json_to_mdl, mtk_simulate, soln_to_csv))
    sim.set_inputs!(args...)
    out = sim.get_outputs!()[1]
    sim.stop!()
    out
end
