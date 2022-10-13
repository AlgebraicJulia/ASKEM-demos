using Catlab
using Catlab.Programs
using Catlab.Present
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.Syntax

@present P(FreeSymmetricMonoidalCategory) begin
    (Petri, Simplex, Trajectory, Number)::Ob
    P::Hom(munit(), Petri)
    multinomial::Hom(Simplex, Simplex)
    solve::Hom(Petri⊗Number⊗Number⊗Simplex, Trajectory)
    measure::Hom(Trajectory, Number)
    sum::Hom(Number⊗Number, Number)
    ptwisesum::Hom(Simplex⊗Simplex, Simplex)
    difference::Hom(Number⊗Number, Number)
end


prog = @program P (θ::Simplex, p::Simplex, t₀::Number, t₁::Number) begin
    θ₁ = multinomial(θ) 
    θ₂ = ptwisesum(θ₁, multinomial(p))
    traj₁ = solve(P(), t₀, t₁, θ₁)
    traj₂ = solve(P(), t₀, t₁, θ₂)
    y₁ = measure(traj₁)
    y₂ = measure(traj₂)
    return difference(y₁, y₂)
end

using JSON

JSON.print(generate_json_acset(prog.diagram), 2)

generate_json_acset(prog.diagram)
