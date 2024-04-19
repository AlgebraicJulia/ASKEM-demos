# Climate models predict a trend of climate instability, including both higher average temperatures along with greater seasonal variations. One of the implications of increasing temperatures is decreasing air density, which affects flight mechanics. Using both historical and predicted data, this script models the increasing distance aircraft require to achieve takeoff on any given 50 km (?) plot on the Earth for a given data. 

# Specifically, this script ingests NetCDF data for lon-lat-temperature data, interpolates it into a spherical mesh representing the Earth, and constructs a decapode to solve an algebraic equation for each point on this mesh. 

# We import the relevant AlgebriacJulia and miscellaneous packages
using ACSets
using CombinatorialSpaces
using DiagrammaticEquations
using DiagrammaticEquations.Deca 
using Decapodes

using Test
using GLMakie
using DataFrames
using ComponentArrays
using GeometryBasics
Point3D = GeometryBasics.Point3{Float64}
using NetCDF
using Plots
using Interpolations
using UnicodePlots
import UnicodePlots: heatmap

export HISTORICAL, SSP370,
historical, ssp370,
hist_lon_bnds, hist_lat_bnds,
size_historical, size_ssp370,
loc_to_idx,
#pode
takeoff_distance

# Runways vary in composition. Suggesting an enhancement for future work, we provide friction coefficients for each runway. `FrictionCoefficients` will be used only when we specify the initial conditions for the decapode, making the simplifying assumption that every point on the sphere has the same friction coefficient of concrete asphalt.
FrictionCoefficients=Dict(
  :ConcreteAsphalt => 0.035,
  :TurfHard        => 0.045,
  :TurfShoftGrass  => 0.0500,
  :TurfHardGrass   => 0.85,
  :GroundSoft      => 0.015
)

# We define a struct for storing location information.
struct Location
  lon::Float32 # for retrieving climate data
  lat::Float32
  material::Symbol
end

# We define a struct for airplane payload.
struct Payload
  type::Symbol
  weight::Float32
end

snakes = Payload(:snakes,30)

# As a reference, here are the lon-lat coordinates for Guam and Little Rock.
guam = Location(144.7937, 13.444304, :ConcreteAsphalt)
littlerock = Location(-92.2895, 34.7464, :GroundSoft)

# We define a struct for our airplane.
mutable struct Aircraft
  weight::Float64 # weight (N)
  max_payload::Float64
  wing_area::Float32 # wing reference area (m²)
  total_thrust::Int # (N)
end

lbs_to_newtons(lbs::Float64) = 4.44822162*lbs

# We define the Boeing 747-8F and Lockheed C-130J
boeing = Aircraft(lbs_to_newtons(470000.0)
  ,lbs_to_newtons(292400.0) ,554.0, 66000)

lockheed = Aircraft(lbs_to_newtons(75560.0)
  ,lbs_to_newtons(42000.0), 162.1
  ,13188 # https://www.military.cz/usa/air/in_service/aircraft/c130/c130_en.htm
  )

planes = [boeing, lockheed]

# We extract and read the historical data
HISTORICAL = "../assets/tasmax_Amon_CanESM5-1_historical_r1i1p1f1_gn_185001-201412.nc"
SSP370 = "../assets/tasmax_Amon_CanESM5_ssp370_r1i1p1f1_gn_201501-210012.nc"

historical = ncread(HISTORICAL, "tasmax")
size_historical = size(historical)

# Likewise for the predicted data
ssp370 = ncread(SSP370, "tasmax")
size_ssp370 = size(ssp370)

dataset = cat(historical, ssp370; dims=3)
dims = size(dataset)

# Construct a decapode based on the takeoff distance formula. Since this form is algebraic, solving it simple requires evaluating it at one time-step. 
takeoff_distance = @decapode begin 
  μ::Constant  # FrictionCoefficient
  τ::Constant  # plane density
  T₀::Constant # static thrust, squared
  Vₛ::Constant # stall speed, squared
  m::Constant  # (g/moles), molar mass
  P::Constant  # Pa, pressure
  V::Constant  # m³, volume (or specific gas constant)
  ρ::Constant  # air density

  T::Form0     # K, temperature
  D::Form0     # takeoff distance

  ∂ₜ(D) == T₀ ./ (2 .* μ .* τ .* (P .* m) ./ (V .* T) .* Vₛ)
end

# Here we define the axes for the visualization
months = 1:dims[3] 
xs = range(0, 180; length=dims[1])
ys = range(-90, 90; length=dims[2])

# The decapode is evaluated on the sphere
sphere = Icosphere(4) |> loadmesh;
sim=evalsim(takeoff_distance)
f=sim(sphere, default_dec_generate)

# We make some assumptions and define a function to calculate plane density.
const molar_mass            = 2.896
const mean_atmos_pa         = 101325
const specific_gas_constant = 287.050

density(plane::Aircraft,load::Float64)=(plane.weight+load*plane.max_payload)/(plane.wing_area) 

# We define a function `meshpolate` which calculates the takeoff distance for a given month and plane. By default, the plane is assumed to be unloaded, but we specify 80% and 100% load.
function meshpolate(data::Matrix{Float32}, plane, load=0.0)
  out = []
  plane ∈ planes || err("Invalid plane")
  linear_interp = linear_interpolation((xs,ys), data) # interpolate
  T = map(sphere[:point]) do pt
      spt = SpherePoint(CartesianPoint(pt))
      linear_interp(
                    isnan(phi(spt)) ? 0 : phi(spt)/pi * 180 + 90, # NAN issue
                    theta(spt) / pi * 180)
  end
  params=(μ=FrictionCoefficients[:ConcreteAsphalt] # 
    ,τ=density(plane,load)    #
    ,T₀=plane.total_thrust    # static thrust, squared
    ,Vₛ=1.0                   # stall speed, squared
    ,m=molar_mass             # molar mass (g/mol)
    ,P=mean_atmos_pa          # pressure (Pa) # TODO this is mean atmospheric pressure
    ,ρ=1.0                    # ρ
    ,V=specific_gas_constant) # volume
  forms=ComponentArray(D=zeros(nv(sphere)), T=T)
  copied=copy(forms)
  f(copied,forms,params,0)
  result = copied.D
end   

# We initialize a dataframe for storing our results
plane_symb = [:boeing, :lockheed]
plane_dict = Dict(:boeing => boeing, :lockheed => lockheed)
df = DataFrame(takeoff_m = Float64[], month = Int[], plane = Symbol[], load = Float64[])

# For each month and plane, calculate and store the takeoff distance assuming 80% load.
for month ∈ months, plane ∈ plane_symb, load ∈ [0.8, 1]
  local out = meshpolate(dataset[:,:,month], plane_dict[plane], load)
  for x in out
    push!(df, [x month plane load])
  end
end

# Here we produce the output.
using Query

fig=Figure()
ax = LScene(fig[1,1], scenekw=(lights=[],))
cam3d!(ax)

sl_load  = Slider(fig[4, 1], range=[0.8, 1], startvalue=0.8)
sl_plane = Slider(fig[3, 1], range=plane_symb, startvalue=:boeing)
sl_time  = Slider(fig[2, 1], range=months, startvalue=1)
sphere_at_time = lift(sl_time.value, sl_plane.value, sl_load.value) do month, plane, load
  @from i in df begin
    @where i.plane == plane && i.month == month && i.load == load
    @select i.takeoff_m
    @collect  
  end
end

msh = mesh!(ax, sphere, color=sphere_at_time, colormap=:jet)
Colorbar(fig[1,2], msh)

display(fig)
