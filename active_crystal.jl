# This Julia code can be use to reproduce the results of the Paper "Lattice-dependent orientational order in active crystals" by Till Welker and Ricard Alert

# The code is structured as follows 

# Definitions of functions. These functions do not need to be adjusted to reproduce results)
#  1. Import Packages
#  2. Functions for Simulations in one and two dimensions
#  3. Read in saved data
#  4. Analysis of order
#  5. Plot Configurations

# Setting Parameters and calling functions. Here you can tune the parameters

#  6. Set Parameters
#  7. Call Functions Here (Example given)




# ---- 1. Import Packages and general functions ----

using CSV
using DataFrames
using Statistics
using PyPlot


function setvis(xlabel = "", ylabel = "",size= (3,2)) # This sets some visual presets for plots
    plt.clf()
    fig, ax = plt.subplots(figsize = (size[1],size[2]))
    fig = plt.gcf()
    fig.patch.set_alpha(0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="y", direction="in", which = "both")
    ax.tick_params(axis="x", direction="in", which = "both")
    return fig,ax
end

# ---- 2. Functions for Simulations in one and two dimensions ----

function get_torque1D(θ,Γ0,Ω) # This function calculates the torques. It is called by the simulation1D method
    θ_l = circshift(θ, 1) # List of left neighbours
    θ_r = circshift(θ, -1) # List of right neighbours

    Γ_XY = Γ0*(Ω+1)/2*(sin.(θ_r.-θ) .+ sin.(θ_l.-θ)) # XY alingment 
    Γ_LA = -Γ0*(Ω-1) *sin.(-2 .* θ) # Lattice alignment
    Γ_MA = Γ0*(Ω-1)/2*(sin.(-θ .-θ_r) .+ sin.(-θ .-θ_l)) # Mirror alignment
    return Γ_XY+Γ_LA+Γ_MA
end

function simulation1D(Γ0,Ω,n,dt,ttot,tsave) # Run for simulation in one dimension. Set the parameters below and call the function to run simulation
    θ = 2*pi*rand(n)    # Random initial condition
    mkdir("1D_Gamma$(Γ0)_Omega$(Ω)")
    for it in 1:ttot    # Iteration in time
        torque = get_torque1D(θ,Γ0,Ω)    # Calucalte Torques
        θ = θ .+ torque .*dt .+ sqrt(2*dt).* randn(n) # Updates Position
        if it%tsave == 0    # Save Data
            counter = round(Int, it/tsave)
            CSV.write("1D_Gamma$(Γ0)_Omega$(Ω)/$counter.csv",  Tables.table(θ), writeheader=false)
            print("$(it/ttot*100)% ")
        end
    end
end


function get_torque2D_square(θ,Γ0,Ω)   # This function calculates the torques. It is called by the simulation2D method
    θ_l = circshift(θ, (0, 1)) # left neighbours
    θ_r = circshift(θ, (0, -1)) # right neighbours
    θ_d = circshift(θ, (-1, 0)) # downwards neighbours
    θ_u = circshift(θ, (1, 0)) # upwards neighbours

    Γ_XY = Γ0*(Ω+1)/2*(sin.(θ_l.-θ).+sin.(θ_r.-θ).+sin.(θ_d.-θ).+sin.(θ_u.-θ)) # XY alignment
    Γ_MA = Γ0*(Ω-1)/2*(sin.(-θ_l.-θ).+sin.(-θ_r.-θ).-sin.(-θ_d.-θ).-sin.(-θ_u.-θ)) # Mirror alignment on a square lattice
    return Γ_XY .+ Γ_MA
end



function simulation2D_square(Γ0,Ω,n,dt,ttot,tsave) # Run for simulation in two dimension. Set the parameters below and call the function to run simulation
    θ = 2*pi*rand(n,n)  # Random initial condition
    mkdir("2D_square_Gamma$(Γ0)_Omega$(Ω)")
    for it in 1:ttot  # Iteration in time
        Γ = get_torque2D_square(θ,Γ0,Ω) # Calucalte Torques
        θ = θ .+ Γ .*dt .+ sqrt(2*dt).* randn(n) # Updates Position
        if it%tsave == 0 # Save Data
            counter = round(Int, it/tsave)
            CSV.write("2D_square_Gamma$(Γ0)_Omega$(Ω)/$counter.csv",  Tables.table(θ), writeheader=false)
            print("$(it/ttot*100)% ")
        end
    end
end


function get_torque2D_triangular(θ,Γ0,Ω)   # This function calculates the torques. It is called by the simulation2D method
    # the lattice sides are indexed using the “double-width” convention given in https://www.redblobgames.com/grids/hexagons/ note that here the first index refers to the line and the second to the row

    θ_l = circshift(θ, (0, 2)) # left neighbours
    θ_r = circshift(θ, (0, -2)) # right neighbours
    θ_lu = circshift(θ, (-1, 1)) # left-upwards neighbours
    θ_ld = circshift(θ, (1, 1)) # left-downwards neighbours
    θ_ru = circshift(θ, (-1, -1)) # right-upwards neighbours
    θ_rd = circshift(θ, (1, -1)) # right-downwards neighbours

    Γ_XY = Γ0*(Ω+1)/2*(sin.(θ_l.-θ).+sin.(θ_r.-θ).+sin.(θ_ld.-θ).+sin.(θ_lu.-θ).+sin.(θ_rd.-θ).+sin.(θ_ru.-θ)) # XY alignment
    Γ_MA = Γ0*(Ω-1)/2*(sin.(-θ_l.-θ).+sin.(-θ_r.-θ).+sin.(2*pi/3 .- θ_ld.-θ).+sin.(4*pi/3 .- θ_lu.-θ).+sin.(4*pi/3 .-θ_rd.-θ).+sin.(2*pi/3 .- θ_ru.-θ)) # Mirror alignment on a triangular lattice

    return Γ_XY .+ Γ_MA
end



function simulation2D_triangular(Γ0,Ω,n,dt,ttot,tsave) # Run for simulation in two dimension. Set the parameters below and call the function to run simulation
    θ = 2*pi*rand(n,2*n)  # Random initial condition
    mkdir("2D_triangular_Gamma$(Γ0)_Omega$(Ω)")
    for it in 1:ttot  # Iteration in time
        Γ = get_torque2D_triangular(θ,Γ0,Ω) # Calucalte Torques
        θ = θ .+ Γ .*dt .+ sqrt(2*dt).* randn(n) # Updates Position
        if it%tsave == 0 # Save Data
            counter = round(Int, it/tsave)
            CSV.write("2D_triangular_Gamma$(Γ0)_Omega$(Ω)/$counter.csv",  Tables.table(θ), writeheader=false)
            print("$(it/ttot*100)% ")
        end
    end
end


# ---- 3. Read in saved data ----

function readin1D(Γ0,Ω,t) # Read in the "t"th saved file for the simulation with Γ0 and Ω for one dimension
    data = CSV.File("1D_Gamma$(Γ0)_Omega$(Ω)/$t.csv"; header=false) |> DataFrame
    θ = Matrix(data)
    return θ
end

function readin2D_square(Γ0,Ω,t) # Read in the "t"th saved file for the simulation with Γ0 and Ω for two dimension
    data = CSV.File("2D_square_Gamma$(Γ0)_Omega$(Ω)/$t.csv"; header=false) |> DataFrame
    θ = Matrix(data)
    return θ
end

function readin2D_triangular(Γ0,Ω,t) # Read in the "t"th saved file for the simulation with Γ0 and Ω for two dimension
    data = CSV.File("2D_triangular_Gamma$(Γ0)_Omega$(Ω)/$t.csv"; header=false) |> DataFrame
    θ = Matrix(data)
    return θ
end

# ---- 4. Analysis of order ----

function nematic_order(θ) # Calucalte the nematic order of the "t"th saved file for the simulation with Γ0 and Ω
    S = sqrt(mean(cos.(2*θ))^2+mean(sin.(2*θ))^2)
    return S
end

function polar_order(θ) # Calucalte the polar order of the "t"th saved file for the simulation with Γ0 and Ω
    P = sqrt(mean(cos.(θ))^2+mean(sin.(θ))^2)
    return P
end



# ---- 5. Plot Configurations ----

function spin_configuration_1D(Γ0,Ω,t) # Plot configurations of Spins on a chain
    θ = readin1D(Γ0,Ω,t)    # Read in file

    fig,ax = setvis("space","",(5,0.5))     # Set up plot
    plt.xticks([])
    plt.yticks([])
    plt.axis("equal")
    plt.xlim(-2,52)
    plt.ylim(-2,2)

    for i in 1:25    # Iterate through 25 spins
        plt.scatter([2*i],[0],color = "black",s = 2,zorder = 4)
        plt.plot([2*i,2*i+cos(θ[i])],[0,sin(θ[i])],color = plt.get_cmap("hsv")(mod(θ[i],2*pi)/2/pi))  # Plot orientations
    end
    

    plt.savefig("1D_spin_configuration_Gamma$(Γ0)_Omega$(Ω).pdf", dpi = 300,bbox_inches = "tight",pad_inches=0.01) # Save figure
end



function spin_configuration_2D_square(Γ0,Ω,t) # Plot configurations of Spins on a square lattice
    θ = readin2D_square(Γ0,Ω,t)    # Read in file

    setvis("","",(2,2))  # Set up plot
    plt.xlim(0,11)
    plt.ylim(0,11)
    plt.xticks([])
    plt.yticks([])

    for i in 1:10   # Iterate over 10x10 spins
        for j in 1:10
            plt.plot([i,i+0.4*cos(θ[j,i])],[j,j+0.4*sin(θ[j,i])],color = plt.get_cmap("hsv")(mod(θ[i],2*pi)/2/pi)) # Plot orientations
            plt.scatter([i],[j], c="black",s = 5)
        end
    end
    plt.savefig("2D_square_spin_configuration_Gamma$(Γ0)_Omega$(Ω).pdf", dpi = 300,bbox_inches = "tight",pad_inches=0.01) # Save figure
end


function spin_configuration_2D_triangular(Γ0,Ω,t) # Plot configurations of Spins on a triangular lattice
    θ = readin2D_triangular(Γ0,Ω,t)    # Read in file

    setvis("","",(2,2))  # Set up plot and create mash
    plt.xlim(0.5,9.5)
    plt.ylim(sqrt(3)/2*0.5,sqrt(3)/2*9.5)
    plt.xticks([])
    plt.yticks([])
    for i in 1:30
        plt.plot([i,-20+i],[-sqrt(3)*20,0],color = "black",alpha = 0.1)
        plt.plot([-10+i,i],[-sqrt(3)*10,0],color = "black",alpha = 0.1)
        plt.plot([0,12],[-sqrt(3)/2*i,-sqrt(3)/2*i],color = "black",alpha = 0.1)
    end

    for i in 1:20   # Iterate over 10x10 spins. We again use the “double-width” convention given in https://www.redblobgames.com/grids/hexagons/ note that here the first index refers to the line and the second to the row
        for j in 1:10
            if (i+j)%2 != 0
                continue
            end
            if j%2 == 0
                plt.plot([i/2,i/2+0.4*cos(θ[j,i])],[sqrt(3)/2*j,sqrt(3)/2*j+0.4*sin(θ[j,i])],color = plt.get_cmap("hsv")(mod(θ[j,i],2*pi)/2/pi)) # Plot orientations
                plt.scatter([i/2],[sqrt(3)/2*j], c="black",s = 5)
            else
                plt.plot([i/2,i/2+0.4*cos(θ[j,i])],[sqrt(3)/2*j,sqrt(3)/2*j+0.4*sin(θ[j,i])],color = plt.get_cmap("hsv")(mod(θ[j,i],2*pi)/2/pi)) # Plot orientations
                plt.scatter([i/2],[sqrt(3)/2*j], c="black",s = 5)
            end
        end
    end
    plt.savefig("2D_triangular_spin_configuration_Gamma$(Γ0)_Omega$(Ω).pdf", dpi = 300,bbox_inches = "tight",pad_inches=0.01) # Save figure
end

# ---- 6. Set Parameters ----

Γ0 = -10 # This is the rescaled torque "\tilde Γ0 * l/a" in the paper
Ω = -1  # This is the distance dependence parameter

n = 1000 # In 1D: number of spins. In 2D: nxn spins
dt = 0.0001  # Timestep 
ttot = Int(100/dt) # Number of simulated timesteps 
tsave = Int(1/dt) # A file is saved every tsave timesteps


# ---- 7. Call Functions Here (The functions called here are an example) ----

simulation1D(Γ0,Ω,n,dt,ttot,tsave)
spin_configuration_1D(Γ0,Ω,100)
