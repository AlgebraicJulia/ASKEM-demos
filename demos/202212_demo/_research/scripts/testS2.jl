

types = formInfType()

SIRD = formSIRD()
states_sird_aug = sirdAugStates()
SIRD_aug = augLabelledPetriNet(SIRD,states_sird_aug)
SIRD_typed = typeSIRD(SIRD_aug, types)
AlgebraicPetri.Graph(dom(SIRD_typed))

Vax = formVax()
states_vax_aug = vaxAugStates()
Vax_aug = augLabelledPetriNet(Vax,states_vax_aug)
Vax_typed = typeVax(Vax_aug, types)
AlgebraicPetri.Graph(dom(Vax_typed))

MA_aug = makeMultiAge(2)
MA_typed = typeAge(MA_aug,types)
AlgebraicPetri.Graph(dom(MA_typed))

SIRD_MA = typed_stratify(SIRD_typed,MA_typed)
AlgebraicPetri.Graph(dom(SIRD_MA))
SIRD_MA_Vax = typed_stratify(SIRD_MA,Vax_typed)
AlgebraicPetri.Graph(dom(SIRD_MA_Vax))

SIRD_Vax = typed_stratify(SIRD_typed,Vax_typed)
AlgebraicPetri.Graph(dom(SIRD_Vax))
SIRD_Vax_MA = typed_stratify(SIRD_Vax,MA_typed)
AlgebraicPetri.Graph(dom(SIRD_Vax_MA))


incident(dom(sviivr_age2),(:id,:id2),:tname)
test = loadSVIIvR("../data/CHIME_SVIIvR_dynamics_BiLayer.json")
testT = [t_interact, t_interact, t_disease, t_interact, t_interact, t_disease, t_disease]
testI, testO = formInfTypeIandO(test,testT)
test_typed = ACSetTransformation(test, types,
           S = [s, s, s, s, s],
           T = testT,
           I = testI,
           O = testO,
           Name = name -> nothing 
           )
test_age = typed_stratify(test_typed,mage_typ)