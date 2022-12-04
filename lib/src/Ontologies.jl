module Ontologies
export infectious_ontology, vector_borne_ontology, vaccine_manufacture_ontology, strip_names

using AlgebraicPetri

function strip_names(p::LabelledPetriNet)
  map(p, Name = name -> nothing)
end

const infectious_ontology = LabelledPetriNet(
    [:Pop],
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop)
  )

const vector_borne_ontology = LabelledPetriNet(
    [:Host, :Vector],
    :infect=>((:Host, :Vector)=>(:Host, :Vector)),
    :host_disease=>(:Host=>:Host),
    :host_strata=>(:Host=>:Host),
    :vector_disease=>(:Vector=>:Vector),
    :vector_strata=>(:Vector=>:Vector)
  )

const vaccine_manufacture_ontology = LabelledPetriNet(
    [:Pop, :Vaccine],
    :produce=>(() => :Vaccine),
    :vaccinate=>((:Pop,:Vaccine) => :Pop),
    :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    :disease=>(:Pop=>:Pop),
    :strata=>(:Pop=>:Pop)
  )

end
