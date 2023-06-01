### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 915f4df0-bc3a-11ec-2560-c9772e143679
begin 
	using Pkg#, Revise
	Pkg.activate(Base.current_project())
	# Pkg.instantiate()
	using Catlab.CategoricalAlgebra
	using Catlab.Present, Catlab.Theories
	using AlgebraicPetri
	using AlgebraicPetri: Graph
	using Semagrams
end;

# ╔═╡ aad0407a-4979-4ea4-94d3-38f1bb6d3de1
begin 
	using Random
	# rng = Random.default_rng()
	Random.seed!(1234)
end;

# ╔═╡ 60d7decb-a7fa-494e-93ca-26d8b957cfcf
include("lib_stratify.jl");

# ╔═╡ 287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
md"""## Model Stratification"""

# ╔═╡ 4330976d-ceeb-4a75-97d0-a8375ede795b
md"""### Load Disease Model from JSON"""

# ╔═╡ 2cf45122-0246-4a07-9efb-e12c76fa9d98
md"""**Load Disease Model from Original MIRA JSON**"""

# ╔═╡ 9efe6732-77d5-4659-8ccd-4b01583454a3
md"""To read in the MIRA model from the JSON file we need to define a new ACSet Type."""

# ╔═╡ 9806cc7e-f14b-453f-af9a-0c0c87e6559c
begin
	@present TheoryOrigMIRANet <: TheoryLabelledReactionNet begin
    	MID::AttrType
    	MCTX::AttrType
    	Template::AttrType
    	mira_ids::Attr(S, MID)
    	mira_context::Attr(S, MCTX)
    	template_type::Attr(T, Template)
    	parameter_name::Attr(T, Name)
    	parameter_value::Attr(T, Rate)
	end
	@abstract_acset_type AbstractOrigMIRANet <: AbstractLabelledReactionNet
	@acset_type OrigMIRANet(TheoryOrigMIRANet) <: AbstractOrigMIRANet
end;

# ╔═╡ f4dc05d6-73ad-4496-97fc-f498c89b8eb9
begin
	load_mira(fname) = begin
	mdl_orig_mira = read_json_acset(OrigMIRANet{Any,Any,Any,Any,Any,Any}, fname)
	mdl_orig_mira[:sname] .= Symbol.(mdl_orig_mira[:sname])
	mdl_orig_mira[:tname] .= Symbol.(mdl_orig_mira[:tname])
	mdl_orig_mira
end
fname = "BIOMD0000000971_petri_orig.json"
mdl_orig = load_mira(fname)
end

# ╔═╡ 02a79d90-34dc-4d6c-b2c8-388acfc5ab95
md"""The model can then be read in as that type and converted to a LabelledPetriNet."""

# ╔═╡ 818b653e-b750-4bdc-8689-ff15c43b8246
# mdl_orig = load_mira("BIOMD0000000971_petri_orig.json");

# ╔═╡ d3043b81-e96b-4b4a-b13e-ec6b3c09868a
md"""The graph of the selected MIRA disease model in this original format is"""

# ╔═╡ 455cbff2-c001-4427-89ee-87393d978d06
AlgebraicPetri.Graph(LabelledPetriNet(mdl_orig))

# ╔═╡ f0f007c1-6528-4a85-a007-b31881297f6a
md"""**Load Disease Model from Curated MIRA JSON**"""

# ╔═╡ 0f12d42f-1a6a-4ca0-bd60-3cb9f38fa518
md"""The original MIRA format is difficult for visualization, so we curated the JSON file. 
This curated format also required its own ACSet Type."""

# ╔═╡ c7a6ed94-e97a-4a96-8f07-1370e71038b9
begin
	@present TheoryMIRANet <: TheoryLabelledReactionNet begin
    	MID::AttrType
    	MCTX::AttrType
    	Template::AttrType
    	mira_ids::Attr(S, MID)
    	mira_context::Attr(S, MCTX)
    	template_type::Attr(T, Template)
    	parameter_name::Attr(T, Name)
    	# parameter_value::Attr(T, Rate)
	end
	@abstract_acset_type AbstractMIRANet <: AbstractLabelledReactionNet
	@acset_type MIRANet(TheoryMIRANet) <: AbstractMIRANet
end;

# ╔═╡ a523d826-3d49-41eb-af52-849d24f0d919
md"""After reading in, rates and concentrations missing from the model still need to be filled. 
We provide default values."""

# ╔═╡ e2c11399-7878-4602-a276-e190857b3fa6
mdl_disease_mira = load_mira_curated("BIOMD0000000971_petri_curated.json", 0.1);

# ╔═╡ a2369b71-43a6-42bf-b322-82e56da8c151
md"""To perform the stratifications and other computations, we extract the model data into the LabelledReactionNet and LabelledPetriNet types, removing the extra MIRA fields."""

# ╔═╡ 6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
begin
	mdl_disease_rxnet = LabelledReactionNet{Any,Any}()
	copy_parts!(mdl_disease_rxnet,mdl_disease_mira)
end;

# ╔═╡ 560fda7a-4f4e-49e9-bf09-d82557286a83
mdl_disease = LabelledPetriNet(mdl_disease_rxnet);

# ╔═╡ 980ebb40-9fe5-46ac-9a74-d324de2d1b26
md"""The graph of the selected MIRA disease model in the curated format is"""

# ╔═╡ 0e104aa0-07c8-4870-a976-7fc6cf8c25de
AlgebraicPetri.Graph(mdl_disease)

# ╔═╡ 93faa61f-55e6-4ed3-a0dd-7f3b153de71d
ns(mdl_disease), nt(mdl_disease)

# ╔═╡ 2c24347b-9c55-4e14-876a-ea8e2867eaa2
md"""### Define Type System"""

# ╔═╡ 23f214d0-ab7c-4325-8ea5-25591f5d571c
md"""Performing stratification requires specifying the types of the components of the models to be stratified. This requires specification of a type system."""

# ╔═╡ b65b454d-8347-40e0-a4fd-44a353f3fe59
md"""We define a type system with one state and three transitions."""

# ╔═╡ f0767520-04f6-4e97-a889-cf5e45d70b4c
begin
	types′ = LabelledPetriNet([:Pop],
    	:infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
    	:disease=>(:Pop=>:Pop),
    	:strata=>(:Pop=>:Pop))
	types = map(types′, Name=name->nothing)
end;

# ╔═╡ 7c003649-03c0-4793-a673-d380e9d199ae
AlgebraicPetri.Graph(types′)

# ╔═╡ fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
md"""**Define Typed Disease Model**"""

# ╔═╡ 8ddeb97d-eb39-4a9d-9eaf-fa955c9019e5
md"""We can now specify the typing of the MIRA disease model by specifying a homomorphism into the type system model."""

# ╔═╡ 4f0fc639-ae18-4c4f-9bd8-aa90037cb759
md"""In this case, the six interaction transitions are of infect type and all other transitions are of disease type."""

# ╔═╡ a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
begin 
	mdl_disease_typed = homomorphism(mdl_disease, types;
    	initial=(T=[1,1,2,2,2,2,2,2,2,2,2,1,2,1,1,1],),
    	type_components=(Name=x->nothing,))
	@assert is_natural(mdl_disease_typed)
end

# ╔═╡ 12a51952-34ab-47a3-a318-c4f58e3f3c12
md"""### Draw Stratification Model"""

# ╔═╡ c89f40f2-3aa9-45c0-9cc8-0abf0a67564f
md"""Here we use `Semagrams.jl` to form the stratification model.""" 

# ╔═╡ 2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
@bind mdl_strat_sema Semagram{StratPetri}(
	"https://semagrams-builds.s3.amazonaws.com/00b3527/petri/main.js",
	"Petri",
	strat_petri_decoders
)

# ╔═╡ cdcbeb7c-20d3-4fd9-9140-aa2db8649c64
mdl_strat_sema

# ╔═╡ b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
md"""### Compute Stratified Model"""

# ╔═╡ 1185894b-063f-402f-96e2-8852430ac4d5
md"""To compute the stratified model, we must fully specify the available transitions and their types."""

# ╔═╡ 6b4584c6-dafa-466f-9879-299556a51277
md"""**Stratification Structs**"""

# ╔═╡ 9baf00c0-1b37-43f8-8a96-f75f45b33a77
md"""Here, we use a high-level, simplified form to make this specification. In this form, the desired transitions are indicated via stratification structs."""

# ╔═╡ 89199dc4-128b-49be-988c-bc69f5b44d8c
md"""Specifically, for each state of each model, we specify the transition types of the other model in which each state participates."""

# ╔═╡ 74e034c5-8f41-47e0-aea1-bdbe401fecd6
md"""In this case, the SQ, H, EQ, and D states cannot travel between cities, but all other states can. These transitions are of type strata. Both cities participate in disease and infect transitions."""

# ╔═╡ 454f4cb1-1535-4aba-a824-181489d324de
begin
	disease_ss = StrataSpec(mdl_disease_typed, 
		[[:strata],[:strata],[:strata],[:strata],[],[],[:strata],[],[]])
	strat_ss = tostrataspec(mdl_strat_sema, types)
end;

# ╔═╡ 1a35abc0-66cf-41ae-b313-7dd248cef669
md"""**Stratified Model**"""

# ╔═╡ 107561be-dcff-4613-bb78-54dee05ed625
md"""From the stratification structs and the type system, we can form the stratified model."""

# ╔═╡ 2f2d6daf-c0b8-4501-86a5-65a844c9359e
mdl_stratified, obs_stratified = stratify(disease_ss, strat_ss, types′);

# ╔═╡ 7836fa75-5ee2-410b-aa12-eded84fac180
md"""The resulting stratified model is"""

# ╔═╡ dad0382a-3656-4bc9-acf6-b914916b15ac
AlgebraicPetri.Graph(mdl_stratified)

# ╔═╡ cf09b8ec-ee28-4823-b895-cdbb339211ab
nparts(mdl_stratified, :S), nparts(mdl_stratified, :T)

# ╔═╡ 1663ff26-27b4-40e7-95fe-20b8ad5333f5
md"""Note that from the stratification we also obtain a particular observation function given by the projection of the pullback onto the disease model, i.e., the observations are the disease states summed across the corresponding stratified states."""

# ╔═╡ 1c42f5e2-4e5c-4616-ab2d-c4d3ac877b9b
md"""### Save Stratified Model"""

# ╔═╡ 8182d1c3-5114-4707-9442-226a6f239e13
md"""We can now save the stratified model to a JSON file."""

# ╔═╡ cbcc0224-dae5-49c3-81eb-b1d37abbd526
write_json_acset(mdl_stratified, "semagrams_stratified_model.json")

# ╔═╡ 4c349031-5d38-4657-9724-b5480a3e1919
write_json_acset(mdl_strat_sema, "semagrams_stratification_model.json")

# ╔═╡ 6b2b41bd-ea99-4917-8eec-d6963243f35a
md"""## Three Stratified Models"""

# ╔═╡ 4fab7265-0d9e-435e-9aa5-972bbf7f5ab8
md"""Here we demonstrate three stratifications of the SIRD disease model and some of the challenges that arise in general stratification."""

# ╔═╡ 1a5c2273-a1d1-44a4-a312-db8e9cd3c7de
md"""**SIRD Disease Model**"""

# ╔═╡ 3f64b4e2-ba9f-4fde-a8bb-c933bd31db0f
md"""The SIRD disease model is given by"""

# ╔═╡ 44328cee-fbf9-42db-afeb-7fc48815952e
SIRD = LabelledPetriNet([:S, :I, :R, :D],
    :inf => ((:S,:I) => (:I,:I)),
    :recover => (:I=>:R),
    :death => (:I=>:D));

# ╔═╡ 3b8deba2-3514-48b6-b7e0-87afa82f3d6b
md"""and visualized as"""

# ╔═╡ 84fec776-b42c-4066-ad45-624b1eb93460
AlgebraicPetri.Graph(SIRD)

# ╔═╡ 1bdbc909-5582-48ff-b2e2-850871542031
md"""As with the previous models, we also must specify the typing of the SIRD disease model."""

# ╔═╡ d482440e-2cce-4f3d-b4cb-208861a71187
begin 
	SIRD_typed = homomorphism(SIRD, types;
    	initial=(T=[1,2,2],I=[1,2,3,3],O=[1,2,3,3]),
    	type_components=(Name=x->nothing,))
	@assert is_natural(SIRD_typed)
end

# ╔═╡ f7d18d80-f2a6-474c-8115-f323e09637d8
md"""**SIRD-Quarantine Model**"""

# ╔═╡ 47bbdd12-2cd5-4e11-8b68-f0bdb0a0cf33
md"""Here we stratify SIRD with a model of Quarantine."""

# ╔═╡ 70bf3253-249d-4def-bbb0-e6f23eb8d717
md"""The Quarantine model is given by"""

# ╔═╡ c4f93cf2-9b30-44b2-9b7a-beee37686f37
begin
	Quarantine = LabelledPetriNet([:Q,:NQ],
    	:quarantine => ((:NQ)=>(:Q)),
    	:unquarantine => ((:Q)=>(:NQ)));
	Quarantine_typed = homomorphism(Quarantine, types;
    	initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(Quarantine_typed)
end

# ╔═╡ 259c785c-aef9-4363-bd2b-996a38797147
md"""and visualized as"""

# ╔═╡ 3c922e81-61af-4789-b2f4-e96dd981810f
AlgebraicPetri.Graph(Quarantine)

# ╔═╡ 178cf6fb-d04f-4c75-8b4a-2647194d9199
md"""To compute the stratification, we must again specify the desired stratification structs. That is, we must specify the transition types of the other model in which each state participates."""

# ╔═╡ fba7130e-7fbc-4ad2-8f86-ab8ea5b6cbca
md"""In this case, the S, I, and R states can undergo the strata transitions, but D state cannot. Meanwhile, both the Q and NQ states can undergo disease transitions, but Q cannot participate in infection interactions."""

# ╔═╡ cc4d55cc-3ad1-47ed-9755-2ef3c616bf57
begin
	SIRD_Q_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	Quarantine_ss = StrataSpec(Quarantine_typed, [[:disease], [:disease,:infect]])
	mdl_SIRD_Q, obsSIRD_Q = stratify(SIRD_Q_ss, Quarantine_ss, types′);
end;

# ╔═╡ 76452de3-13dd-48eb-852c-f7bb4bbc38fc
md"""The stratified model is thus"""

# ╔═╡ 0de945ad-093d-4975-8e48-a0748e2414b1
AlgebraicPetri.Graph(mdl_SIRD_Q)

# ╔═╡ 385a2551-4733-4fb0-9d5d-9763f6757578
md"""**SIRD-Mask Model (v1 - High-Level, Simplified Stratification Form)**"""

# ╔═╡ e3ca90f4-b754-4834-89ce-81b8eb8ee614
md"""Here we stratify SIRD with a model of Masking."""

# ╔═╡ 2fb2e995-14bf-4323-800e-7c0cfadb4a7a
md"""In our model of masking, non-masked individuals can both infect and become infected, while masked individuals can become infected but do not infect others."""

# ╔═╡ 028fc3bd-8577-4a85-81bd-18f3fec75437
md"""This model is given by"""

# ╔═╡ 452ade8f-6140-4531-9919-51e899b15f33
begin
	Mask = LabelledPetriNet([:M,:NM],
    	:mask => ((:NM)=>(:M)),
    	:unmask => ((:M)=>(:NM)),
		:mu_infect => ((:M,:NM)=>(:M,:NM)),
		:uu_infect => ((:NM,:NM)=>(:NM,:NM)));
	Mask_typed = homomorphism(Mask, types;
    	initial=(T=[3,3,1,1],), type_components=(Name=x->nothing,))
	@assert is_natural(Mask_typed)
end

# ╔═╡ e79d741f-702d-49c9-9e1a-6442cfeb7614
AlgebraicPetri.Graph(Mask)

# ╔═╡ 270ed674-350b-4911-86e2-091fbf619f2a
md"""In forming the stratification structs, we again note that the S, I, and R states can undergo the change in masking, but the D state cannot, while both M and NM need to participate in the disease and infection transitions."""

# ╔═╡ c33580a9-1d22-42a5-b527-018c2448b209
begin
	SIRD_M_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	Mask_ss = StrataSpec(Mask_typed, [[:disease,:infect], [:disease,:infect]])
	mdl_SIRD_M, obsSIRD_M = stratify(SIRD_M_ss, Mask_ss, types′);
end;

# ╔═╡ 3b786190-1a6a-4490-a545-c2e23595d09f
md"""Performing the stratification in this form leads to the following model"""

# ╔═╡ a9043c9e-ec1b-41e8-a617-5a1a4d611415
AlgebraicPetri.Graph(mdl_SIRD_M)

# ╔═╡ 803517ca-ad26-4c3c-b003-78d8b9bdf0cc
md"""**SIRD-Mask Model (v2 - Full-Pullback Stratification Form)**"""

# ╔═╡ 5d1149e0-2007-4738-a658-fe1111147c6d
md"""Unfortunately, undesired transitions are apparent in the stratified model, e.g., a mask susceptible and unmasked susceptible spontaneously forming two infected individuals, ((S,M),(S,NM))=>((I,M),(I,NM))."""

# ╔═╡ 47d72cb1-215d-4f10-88df-d577b3fb77a4
md"""The form for specifying the stratification used above is a high-level, simplified version to avoid needing to augment the component models with identity transitions."""

# ╔═╡ 0e50dc55-a18b-4203-9060-9d8cacc94718
md"""However, the asymmetry between the two roles in the infection interaction and the asymmetry in the participation of masked individuals in that process is not captured by the high-level form of specification and produces anomalous transitions in the stratified model."""

# ╔═╡ c34b47d0-e493-4af7-b14d-89196507022a
md"""This issue can be alleviated by using the full-pullback form for specifying the stratification."""

# ╔═╡ dc591063-59cc-41cb-b194-c033b33f0dfe
md"""This more general form requires augmenting the component models with appropriate identity transitions and specifying the typing of all transitions."""

# ╔═╡ fb274f1d-fab3-4d50-ba87-c8f1cd880822
begin
	# These are definitions for ease of use
	s, = parts(types′, :S)
	t_interact, t_disease, t_strata = parts(types′, :T)
	i_interact1, i_interact2, i_disease, i_strata = parts(types′, :I)
	o_interact1, o_interact2, o_disease, o_strata = parts(types′, :O);
end;

# ╔═╡ 13d1d050-e8ca-49f6-8e7d-23b22efba502
md"""For the SIRD model, identity transitions are added for S, I, and R, but not D. These transitions are of strata type, while the original transitions are still of types interact, disease, and disease, respectively. """

# ╔═╡ 75d3e1cd-b951-4626-9625-9ef1c546ef87
begin 
	SIRD_aug = LabelledPetriNet([:S, :I, :R, :D],
	  :inf => ((:S, :I)=>(:I, :I)),
	  :rec => (:I=>:R),
	  :death => (:I=>:D),
	  :id => (:S => :S),
	  :id => (:I => :I),
	  :id => (:R => :R)
	)
	
	SIRD_aug_typed = ACSetTransformation(SIRD_aug, types,
	  S = [s, s, s, s],
	  T = [t_interact, t_disease, t_disease, t_strata, t_strata, t_strata],
	  I = [i_interact1, i_interact2, i_disease, i_disease, i_strata, i_strata, i_strata],
	  O = [o_interact1, o_interact2, o_disease, o_disease, o_strata, o_strata, o_strata],
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
	)
	
	@assert is_natural(SIRD_aug_typed)
end;

# ╔═╡ 047d247e-2ed2-41ef-b6ad-069f2a92d227
md"""The augmented SIRD model is thus"""

# ╔═╡ c0458439-2af8-4e42-90c7-41f562d0750e
AlgebraicPetri.Graph(dom(SIRD_aug_typed)) # Graph_typed

# ╔═╡ a7756694-a062-43cd-913b-b350beca1b59
md"""For the Mask model, identity transitions are added for both M and NM. These transitions are of disease type, while the original transitions are of types interact, interact, strata, and strata, respectively."""

# ╔═╡ 488e5c05-39ee-4635-9346-730bea5aa5f5
begin
	Mask_aug = LabelledPetriNet([:M,:NM],
    	:mu_infect => ((:M,:NM)=>(:M,:NM)),
		:uu_infect => ((:NM,:NM)=>(:NM,:NM)),
	  	:mask => ((:NM)=>(:M)),
    	:unmask => ((:M)=>(:NM)),
		:id => (:M => :M),
	  	:id => (:NM => :NM)
	)
	
	Mask_aug_typed = ACSetTransformation(Mask_aug, types,
	  S = [s, s],
	  T = [t_interact, t_interact, t_strata, t_strata, t_disease, t_disease],
	  I = [i_interact1, i_interact2, i_interact1, i_interact2, i_strata, i_strata, i_disease, i_disease],
	  O = [o_interact1, o_interact2, o_interact1, o_interact2, o_strata, o_strata, o_disease, o_disease],
	  Name = name -> nothing # specify the mapping for the loose ACSet transform
	)

	@assert is_natural(Mask_aug_typed)
end;

# ╔═╡ eb366215-ab33-446a-aef5-d9d95c93e49d
md"""The augmented Mask model is thus"""

# ╔═╡ ece17ed4-dc74-4bfc-b440-6f34cd560cab
AlgebraicPetri.Graph(dom(Mask_aug_typed)) # Graph_typed

# ╔═╡ e64fadee-e120-4eb6-9d1b-c69b710900cc
md"""With the augmented models, the stratification is given by the pullback."""

# ╔═╡ 06636e60-e7dc-43ee-8541-190881085d16
begin
	typed_stratify(typed_model1, typed_model2) =
		Theories.compose(proj1(CategoricalAlgebra.pullback(typed_model1, typed_model2)), typed_model1)
	mdl2_SIRD_M = typed_stratify(SIRD_aug_typed, Mask_aug_typed)
end;

# ╔═╡ 5aa89767-cd1e-49f5-838d-9efd8bd39ef3
md"""The fully specified stratification leads to the following model, which is correct."""

# ╔═╡ 03ad5137-c6f7-424d-99fb-65ac7f8a2ed2
AlgebraicPetri.Graph(dom(mdl2_SIRD_M)) # Graph_typed

# ╔═╡ d9fd32c2-ed42-4af6-90a7-d144296e222d
md"""**SIRD-Two-City Model**"""

# ╔═╡ 44e9d001-e0df-433a-a05a-00c57e4c97c5
md"""Here we stratify SIRD with the model of movement between Two Cities."""

# ╔═╡ 0cca8978-d7b4-496b-8e48-3e5bf14a4028
md"""As with the MIRA disease model, this stratification can be specified using the previous high-level form, with the Two-City model given by"""

# ╔═╡ 4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
begin
	TwoCity = LabelledPetriNet([:City1,:City2],
	    :travel12 => ((:City1)=>(:City2)),
    	:travel21 => ((:City2)=>(:City1)))
	TwoCity_typed = homomorphism(TwoCity, types;
	    initial=(T=[3,3],), type_components=(Name=x->nothing,))
	@assert is_natural(TwoCity_typed)
end

# ╔═╡ 25ba5ac6-9c41-4dfc-93b2-75664ea59565
AlgebraicPetri.Graph(TwoCity)

# ╔═╡ 38f0c179-16b2-4c4d-9975-8fb4c3f9f97b
begin
	SIRD_TC_ss = StrataSpec(SIRD_typed, [[:strata],[:strata],[:strata],[]])
	TwoCity_ss = StrataSpec(TwoCity_typed, [[:disease,:infect], [:disease,:infect]])
	mdl_SIRD_TC, obsSIRD_TC = stratify(SIRD_TC_ss, TwoCity_ss, types′);
end;

# ╔═╡ 91ad223d-c245-4aca-89ab-cd6ee8c6652d
md"""Resulting in the following stratified model"""

# ╔═╡ 015fa0fd-5e9f-469a-bbb8-c388b946b3ab
AlgebraicPetri.Graph(mdl_SIRD_TC)

# ╔═╡ Cell order:
# ╠═915f4df0-bc3a-11ec-2560-c9772e143679
# ╠═60d7decb-a7fa-494e-93ca-26d8b957cfcf
# ╠═aad0407a-4979-4ea4-94d3-38f1bb6d3de1
# ╟─287a44bb-e4d9-4a8b-b383-13d9c01e6e1e
# ╟─4330976d-ceeb-4a75-97d0-a8375ede795b
# ╟─2cf45122-0246-4a07-9efb-e12c76fa9d98
# ╟─9efe6732-77d5-4659-8ccd-4b01583454a3
# ╠═9806cc7e-f14b-453f-af9a-0c0c87e6559c
# ╠═f4dc05d6-73ad-4496-97fc-f498c89b8eb9
# ╟─02a79d90-34dc-4d6c-b2c8-388acfc5ab95
# ╠═818b653e-b750-4bdc-8689-ff15c43b8246
# ╟─d3043b81-e96b-4b4a-b13e-ec6b3c09868a
# ╠═455cbff2-c001-4427-89ee-87393d978d06
# ╟─f0f007c1-6528-4a85-a007-b31881297f6a
# ╟─0f12d42f-1a6a-4ca0-bd60-3cb9f38fa518
# ╠═c7a6ed94-e97a-4a96-8f07-1370e71038b9
# ╟─a523d826-3d49-41eb-af52-849d24f0d919
# ╠═e2c11399-7878-4602-a276-e190857b3fa6
# ╟─a2369b71-43a6-42bf-b322-82e56da8c151
# ╠═6dec59bf-1fa0-4ad6-afdd-26ec2b92230a
# ╠═560fda7a-4f4e-49e9-bf09-d82557286a83
# ╟─980ebb40-9fe5-46ac-9a74-d324de2d1b26
# ╠═0e104aa0-07c8-4870-a976-7fc6cf8c25de
# ╠═93faa61f-55e6-4ed3-a0dd-7f3b153de71d
# ╟─2c24347b-9c55-4e14-876a-ea8e2867eaa2
# ╟─23f214d0-ab7c-4325-8ea5-25591f5d571c
# ╟─b65b454d-8347-40e0-a4fd-44a353f3fe59
# ╠═f0767520-04f6-4e97-a889-cf5e45d70b4c
# ╠═7c003649-03c0-4793-a673-d380e9d199ae
# ╟─fbe5adf1-1d6b-40a6-b36e-49e3fc76057e
# ╟─8ddeb97d-eb39-4a9d-9eaf-fa955c9019e5
# ╟─4f0fc639-ae18-4c4f-9bd8-aa90037cb759
# ╠═a3eca5f4-21c0-4147-b7d0-19ba3fd7a776
# ╟─12a51952-34ab-47a3-a318-c4f58e3f3c12
# ╟─c89f40f2-3aa9-45c0-9cc8-0abf0a67564f
# ╠═2cc7237d-3a8b-4b86-bb08-6cba2e9e9531
# ╠═cdcbeb7c-20d3-4fd9-9140-aa2db8649c64
# ╟─b598e5b8-38e2-4d7d-afc9-2b0e7b1d2be8
# ╟─1185894b-063f-402f-96e2-8852430ac4d5
# ╟─6b4584c6-dafa-466f-9879-299556a51277
# ╟─9baf00c0-1b37-43f8-8a96-f75f45b33a77
# ╟─89199dc4-128b-49be-988c-bc69f5b44d8c
# ╟─74e034c5-8f41-47e0-aea1-bdbe401fecd6
# ╠═454f4cb1-1535-4aba-a824-181489d324de
# ╟─1a35abc0-66cf-41ae-b313-7dd248cef669
# ╟─107561be-dcff-4613-bb78-54dee05ed625
# ╠═2f2d6daf-c0b8-4501-86a5-65a844c9359e
# ╟─7836fa75-5ee2-410b-aa12-eded84fac180
# ╠═dad0382a-3656-4bc9-acf6-b914916b15ac
# ╠═cf09b8ec-ee28-4823-b895-cdbb339211ab
# ╟─1663ff26-27b4-40e7-95fe-20b8ad5333f5
# ╟─1c42f5e2-4e5c-4616-ab2d-c4d3ac877b9b
# ╟─8182d1c3-5114-4707-9442-226a6f239e13
# ╠═cbcc0224-dae5-49c3-81eb-b1d37abbd526
# ╠═4c349031-5d38-4657-9724-b5480a3e1919
# ╟─6b2b41bd-ea99-4917-8eec-d6963243f35a
# ╟─4fab7265-0d9e-435e-9aa5-972bbf7f5ab8
# ╟─1a5c2273-a1d1-44a4-a312-db8e9cd3c7de
# ╟─3f64b4e2-ba9f-4fde-a8bb-c933bd31db0f
# ╠═44328cee-fbf9-42db-afeb-7fc48815952e
# ╟─3b8deba2-3514-48b6-b7e0-87afa82f3d6b
# ╠═84fec776-b42c-4066-ad45-624b1eb93460
# ╟─1bdbc909-5582-48ff-b2e2-850871542031
# ╠═d482440e-2cce-4f3d-b4cb-208861a71187
# ╟─f7d18d80-f2a6-474c-8115-f323e09637d8
# ╟─47bbdd12-2cd5-4e11-8b68-f0bdb0a0cf33
# ╟─70bf3253-249d-4def-bbb0-e6f23eb8d717
# ╠═c4f93cf2-9b30-44b2-9b7a-beee37686f37
# ╟─259c785c-aef9-4363-bd2b-996a38797147
# ╠═3c922e81-61af-4789-b2f4-e96dd981810f
# ╟─178cf6fb-d04f-4c75-8b4a-2647194d9199
# ╟─fba7130e-7fbc-4ad2-8f86-ab8ea5b6cbca
# ╠═cc4d55cc-3ad1-47ed-9755-2ef3c616bf57
# ╟─76452de3-13dd-48eb-852c-f7bb4bbc38fc
# ╠═0de945ad-093d-4975-8e48-a0748e2414b1
# ╟─385a2551-4733-4fb0-9d5d-9763f6757578
# ╟─e3ca90f4-b754-4834-89ce-81b8eb8ee614
# ╠═2fb2e995-14bf-4323-800e-7c0cfadb4a7a
# ╟─028fc3bd-8577-4a85-81bd-18f3fec75437
# ╠═452ade8f-6140-4531-9919-51e899b15f33
# ╠═e79d741f-702d-49c9-9e1a-6442cfeb7614
# ╟─270ed674-350b-4911-86e2-091fbf619f2a
# ╠═c33580a9-1d22-42a5-b527-018c2448b209
# ╟─3b786190-1a6a-4490-a545-c2e23595d09f
# ╠═a9043c9e-ec1b-41e8-a617-5a1a4d611415
# ╟─803517ca-ad26-4c3c-b003-78d8b9bdf0cc
# ╟─5d1149e0-2007-4738-a658-fe1111147c6d
# ╟─47d72cb1-215d-4f10-88df-d577b3fb77a4
# ╟─0e50dc55-a18b-4203-9060-9d8cacc94718
# ╟─c34b47d0-e493-4af7-b14d-89196507022a
# ╟─dc591063-59cc-41cb-b194-c033b33f0dfe
# ╠═fb274f1d-fab3-4d50-ba87-c8f1cd880822
# ╟─13d1d050-e8ca-49f6-8e7d-23b22efba502
# ╠═75d3e1cd-b951-4626-9625-9ef1c546ef87
# ╟─047d247e-2ed2-41ef-b6ad-069f2a92d227
# ╠═c0458439-2af8-4e42-90c7-41f562d0750e
# ╟─a7756694-a062-43cd-913b-b350beca1b59
# ╠═488e5c05-39ee-4635-9346-730bea5aa5f5
# ╟─eb366215-ab33-446a-aef5-d9d95c93e49d
# ╠═ece17ed4-dc74-4bfc-b440-6f34cd560cab
# ╟─e64fadee-e120-4eb6-9d1b-c69b710900cc
# ╠═06636e60-e7dc-43ee-8541-190881085d16
# ╠═5aa89767-cd1e-49f5-838d-9efd8bd39ef3
# ╠═03ad5137-c6f7-424d-99fb-65ac7f8a2ed2
# ╟─d9fd32c2-ed42-4af6-90a7-d144296e222d
# ╟─44e9d001-e0df-433a-a05a-00c57e4c97c5
# ╟─0cca8978-d7b4-496b-8e48-3e5bf14a4028
# ╠═4c7d9df5-8f92-4fe1-a2e3-82f607d28ccb
# ╠═25ba5ac6-9c41-4dfc-93b2-75664ea59565
# ╠═38f0c179-16b2-4c4d-9975-8fb4c3f9f97b
# ╟─91ad223d-c245-4aca-89ab-cd6ee8c6652d
# ╠═015fa0fd-5e9f-469a-bbb8-c388b946b3ab
