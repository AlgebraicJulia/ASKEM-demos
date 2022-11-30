
using AlgebraicPetri: vectorfield

# AlgebraicDynamics: ContinuousResourceSharer

valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end
AlgebraicPetri.vectorfield(pn::AbstractPetriNet) = begin
    tm = TransitionMatrices(pn)
    dt = tm.output - tm.input
    f(du,u,p,t) = begin
      rates = zeros(eltype(du),nt(pn))
      # u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
      # p_m = [p[tname(pn, i)] for i in 1:nt(pn)]
      u_m = u
      p_m = p
      for i in 1:nt(pn)
        rates[i] = valueat(p_m[i],u,t) * prod(u_m[j] ^ tm.input[i,j] for j in 1:ns(pn))
      end
      for j in 1:ns(pn)
        # du[sname(pn, j)] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
        du[j] = sum(rates[i] * dt[i,j] for i in 1:nt(pn); init = 0.0)
      end
      return du
    end
    return f
end

