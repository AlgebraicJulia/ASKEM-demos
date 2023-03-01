# Scenario 3

## 3.1 Base SIR model


```julia
using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab.CategoricalAlgebra, Catlab.WiringDiagrams

const infectious_ontology = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

Graph(infectious_ontology)
```


```julia
using Catlab.Programs, Catlab.Graphics

sir_uwd = @relation (S,I,R) where (S::Pop, I::Pop, R::Pop) begin
  infect(S,I,I,I) # inf
  disease(I,R) # rec
end

to_graphviz(sir_uwd, box_labels=:name, junction_labels=:variable)
```


```julia
base_names = [:inf, :rec]
typed_sir = oapply_typed(infectious_ontology, sir_uwd, base_names)

Graph(dom(typed_sir))
```

## 3.2 Hospitalization and death

### SIRD

SIR model with death


```julia
with_death_uwd = @relation (S,I,R,D) where (S::Pop, I::Pop, R::Pop, D::Pop) begin
  base(S,I,R)
  disease(I,D) # die
end

sird_uwd = ocompose(with_death_uwd, 1, sir_uwd)

to_graphviz(sird_uwd, box_labels=:name, junction_labels=:variable)
```


```julia
typed_sird = oapply_typed(infectious_ontology, sird_uwd, [base_names; :die])

Graph(dom(typed_sird))
```


```julia
write_json_acset(dom(typed_sird), "scenario3_sird.json")
```

### SIRH

SIR model with hospitalization


```julia
with_hospital_uwd = @relation (S,I,R,H) where (S::Pop, I::Pop, R::Pop, H::Pop) begin
  base(S,I,R)
  disease(I,H) # hosp
  disease(H,R) # hosp_rec
end

sirh_uwd = ocompose(with_hospital_uwd, 1, sir_uwd)

to_graphviz(sirh_uwd, box_labels=:name, junction_labels=:variable)
```


```julia
typed_sirh = oapply_typed(infectious_ontology, sirh_uwd,
  [base_names; [:hosp, :hosp_rec]])

Graph(dom(typed_sirh))
```


```julia
write_json_acset(dom(typed_sirh), "scenario3_sirh.json")
```

### SIRHD

SIR model with both hospitalization and death


```julia
with_hospital_death_uwd = @relation (S,I,R,H,D) where (S::Pop, I::Pop, R::Pop, H::Pop, D::Pop) begin
  base(S,I,R,H)
  disease(I,D) # die
  disease(H,D) # hosp_die
end

sirhd_uwd = ocompose(with_hospital_death_uwd, 1, sirh_uwd)

to_graphviz(sirhd_uwd, box_labels=:name, junction_labels=:variable,
  edge_attrs=Dict(:len => "0.75"))
```


```julia
typed_sirhd = oapply_typed(infectious_ontology, sirhd_uwd,
  [base_names; [:hosp, :hosp_rec, :die, :hosp_die]])

Graph(dom(typed_sirhd))
```


```julia
write_json_acset(dom(typed_sirhd), "scenario3_sirhd.json")
```

## 3.3 Vaccination

SIRHD model with vaccination.


```julia
vaccination_uwd = @relation () where (U::Pop, V::Pop) begin
  strata(U,V) # vac
  infect(U,U,U,U)
  infect(U,V,U,V)
  infect(V,U,V,U)
  infect(V,V,V,V)
end

typed_vaccination = oapply_typed(infectious_ontology, vaccination_uwd,
  [:vac, :UU, :UV, :VU, :VV])

to_graphviz(vaccination_uwd, box_labels=:name, junction_labels=:variable,
  edge_attrs=Dict(:len => "1"), graph_attrs=Dict(:start => "2"))
```


```julia
typed_sirhd_aug = add_reflexives(
  typed_sirhd, [[:strata],[],[],[],[]], infectious_ontology)

typed_vaccination_aug = add_reflexives(
  typed_vaccination,
  [[:disease], [:disease]],
  infectious_ontology
)

typed_sirhd_vac = typed_product(typed_sirhd_aug, typed_vaccination_aug)

Graph(dom(typed_sirhd_vac))
```


```julia
write_json_acset(dom(typed_sirhd_vac), "scenario3_sirhd_vac.json")
```

## 3.4 Age stratification

SIRHD model with vaccination *and* age stratification.


```julia
for n in (2, 8)
  names = [Symbol("Age$i") for i in 1:n]
  typed_age = pairwise_id_typed_petri(infectious_ontology, :Pop, :infect, names)

  typed_age_aug = add_reflexives(
    typed_age,
    repeat([[:disease, :strata]], n),
    infectious_ontology
  )

  typed_sirhd_vac_age = typed_product(typed_sirhd_vac, typed_age_aug)
  net = dom(typed_sirhd_vac_age)
  
  write_json_acset(net, "scenario3_sirhd_vac_age$n.json")
  
  open("scenario3_sirhd_vac_age$n.svg", "w") do io
    show(io, MIME("image/svg+xml"), Graph(net))
  end
end
```

## Bonus: Testing of infectives

### SIRT

SIR model with tested and untested infectives


```julia
sirt_uwd = @relation (S,Iᵤ,Iₜ,R) where (S::Pop, Iᵤ::Pop, Iₜ::Pop, R::Pop) begin
  infect(S,Iᵤ,Iᵤ,Iᵤ) # infᵤ
  infect(S,Iₜ,Iᵤ,Iₜ) # infₜ
  disease(Iᵤ,R) # recᵤ
  disease(Iₜ,R) # recₜ
  disease(Iᵤ,Iₜ) # test
end

sirt_names = [:infᵤ, :infₜ, :recᵤ, :recₜ, :test]
typed_sirt = oapply_typed(infectious_ontology, sirt_uwd, sirt_names)

Graph(dom(typed_sirt))
```


```julia
write_json_acset(dom(typed_sirt), "scenario3_sirt.json")
```

### Testing and vaccination

Add vaccination to this model using stratification.


```julia
typed_sirt_aug = add_reflexives(
  typed_sirt, [[:strata],[],[],[]], infectious_ontology)

typed_sirt_vac = typed_product(typed_sirt_aug, typed_vaccination_aug)

Graph(dom(typed_sirt_vac))
```


```julia
write_json_acset(dom(typed_sirt_vac), "scenario3_sirt_vac.json")
```
