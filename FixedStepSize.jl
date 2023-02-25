using LinearAlgebra
using LinearOperators
using Plots
using MatrixDepot
using SparseArrays
#-----------------------------------------------------
struct Problem
    c::Vector{Float64}
    A::Matrix{Float64} 
    b::Vector{Float64}
    Eqin::Vector{Float64}
    function Problem(c::Vector{Float64}, A::Matrix{Float64}  ,b::Vector{Float64} ,Eqin::Vector{Float64})
        new(c,A,b,Eqin)     
    end
end
function preSolve(problem::Problem)
    m,n = size(problem.A)
    minit = m
    ninit = n
    A_temp = deepcopy(problem.A)
    c_temp= deepcopy(problem.c)
    id=[]
    for i= 1:m
        if (problem.Eqin[i]!=0)
            n=n+1
            append!(id,i)
        end
    end
    m,n = size(problem.A)
    for i = 1:length(id)
        if (Eqin[id[i]]==-1) #<= constraint
            ei= zeros(m,1)   
            ei[id[i]]=1
            A_temp= [A_temp ei]
            append!(c_temp, 0)
        else
            ei= zeros(m,1)
            ei[id[i]] = -1
            A_temp = [A_temp ei] 
            append!(c_temp, 0)
            #should equate Aproblem
        end
    end
    print(c_temp)
    Pr = Problem(c_temp, A_temp, problem.b, problem.Eqin)
    
end
print("Enter number of variables: \n")
n = parse(Int,readline(stdin))
print("Enter number of constraints: \n")
m = parse(Int,readline(stdin))
#m=3
#n=2
A = Matrix{Float64}(undef,m,n)
b = Array{Float64}(undef,m)
c = Array{Float64}(undef,n)
Eqin = Array{Float64}(undef,m)
for i= 1:m
    print("--- Next constraint ---")
    for j = 1:n
         println("Coefficient of X", j,": ")
         X_coef = readline()
         X_coef = parse(Float64, X_coef)
         A[i,j]= X_coef
     end
     print("Enter the type of constraint")
     X_coef = readline()
     X_coef = parse(Int64, X_coef)
     Eqin[i]= X_coef
     print("The solution of the constraint: ")
     X_coef = readline()
     X_coef = parse(Float64, X_coef)
     b[i]= X_coef
 end
 print("--- Objective function ---")
 for i= 1:n
     println("Enter Coefficient of X",i,": ")
     X_coef = readline()
     X_coef = parse(Float64, X_coef)
     c[i]= X_coef
 end
 print("Enter 1- for Maximization and 2- for Minimization: ")
 MinMaxLP = readline()
 MinMaxLP = parse(Int, MinMaxLP)
#delete here
#A = [1.0 1.0 1.0]
#c = [-1.1, 1.0, 0.0]
#b = [6.0]
#Eqin = [0.0,1.0,-1.0]
#MinMaxLP=2
if MinMaxLP == 1
    c .*= -1.0
end
println(c) 
Pr = Problem(c,A,b,Eqin)  
P=preSolve(Pr)
print("\nPost Presolving:     ",P.A)
xsol= []
fval=0
exitflag=0
itr=0
maxIter=100
tol= 1e-8
etaMin =0.995
scalingtec= 6
infeasible =1
c0 = 0

m,n = size(P.A)

if infeasible != 1
    print("The LP problem is infeasible \n")
    exitflag=1
    return
end

x = P.A' * (P.A * P.A')^-1 *P.b
print("A: ",P.A)
print("\nb: ",P.b)
print("\nc: ",P.c)
w = (P.A * P.A')^(-1) *P.A *P.c

print("\n\nw:",w)
s = Matrix{Float64}(undef,1,m)
s = P.c - P.A' * w
#print("\ns: ",s)
delta_x = maximum([-1.5 * minimum(x), 0])
delta_s = maximum([-1.5 * minimum(s), 0])
e = ones(n)

delta_x_s = 0.5 * (x + delta_x * e)' * (s + delta_s * e)
delta_x_c = delta_x + delta_x_s /(sum(s) + n * delta_s)
print("\ne ",delta_x_c)
delta_s_c = delta_s + delta_x_s /(sum(x) + n * delta_x)

x = x + delta_x_c * e
print("\n\nx:",x)
s = s + delta_s_c * e
print("\n\ns:",s)
α = 0.9 
σ = 0.5 #belongs to [0,1]
######################
##                 ###
##    the loop     ###
##                 ###
######################
iterl = []
FArray = Float64[]
xArray=Vector{Float64}[]
sArray=Vector{Float64}[]

for itr= 1:maxIter
    
     itr = itr + 1
     push!(iterl,itr)
     print("\n----------------------------- Iteration: ",itr," --------------------------------------------------")
     
     identity = Matrix(1.0I,n,n)
     augmented_subsystem_1 = hcat(zeros(n,n), P.A', identity)
    # print("\n ASs1: ",augmented_subsystem_1)
     augmented_subsystem_2 = hcat(P.A, zeros(m,m), zeros(m,n))
    # print("\n ASs2: ",augmented_subsystem_2)
     augmented_subsystem_3 = hcat(Diagonal(s), zeros(n), Diagonal(x))
    # print("\n ASs3: ",augmented_subsystem_3)
     augmented_system = vcat(augmented_subsystem_1, augmented_subsystem_2, augmented_subsystem_3)
   #  print("\nASS: ",augmented_system)
     rp = P.A * x .- b
    # print("\nrp: ", rp)
     μ  = x' * s /n
    # print("\nμ : ", μ )
     rc = P.A' * w + s - P.c
    # print("\nrc: ", rc)   
     rLast = Diagonal(x) * Diagonal(s)* ones(n).- μ  .*σ
    # print("\nLast: ", rLast)
    
    if max(μ , norm(rp),norm(rLast)) <= tol
        xsol = x
        fval = P.c' * x
        if MinMaxLP == 1  #Maximization
             fval = -(fval + c0)
        else
             fval = fval + c0
        end
        exitflag = 0
        print("\n-----------WE REACHED AN OPTIMAL POINT ! :D-----------")
        break
        
    end
     AB= vcat(rc .*-1,rp .*-1,rLast .*-1)
   #  print("\nAB: ",AB)
     Deltas = (augmented_system)^-1 * AB
   #  print("\n\nDeltas: ", Deltas)
  
     x[1:length(x)] += α*Deltas[1:length(x)]
     print("\n x: ",x)
     push!(xArray,x)
     w[1:length(w)] += α*Deltas[(1:length(w)).+length(x)]
     print("\n w: ",w)
     xw = length(w)+length(x)
     s[1:length(s)] += α*Deltas[(1:length(s)).+xw]
     print("\n s: ",s)
     push!(sArray,s)
     if mod(itr,1)==0
        if MinMaxLP == 1  #Maximization
             fval = -(P.c'* x + c0)
             push!(FArray,fval)
        else
             fval = P.c'* x + c0
             push!(FArray,fval)
        end
           print("\n objective function: ",fval)
    end
end

 print("\n-------- Using the Central-path fixed step-size Method -------")
    print("\nx*: ",xsol)
    print("\nw*: ",w)
    print("\ns*: ",s)
    z = P.c' * x
    push!(FArray, z)
    print("\n Optimal Solution:", z)

    print("\n\n\n x:",xArray)
    print("\n\n\n s:",sArray)
#     print("\n\n\n f:",FArray)
   # print("\n\n\n iterl:",iterl)
    # PlotAll(xArray,sArray,iterl,FArray)  
    SX1 = Float64[]
    SX2 = Float64[]
    zx = Float64[]
    xz = Float64[]
    lsx=map(v->v[1],xArray)
    print("\n lsx: ",lsx)
    lxz=map(v->v[2],xArray)
    print("\n xz: ",xz)
    for i=1:length(xArray)
        zx=xArray[i][1].*sArray[i][1]
        push!(SX1,zx)
        xz=xArray[i][2].*sArray[i][2]
        push!(SX2,xz)
    end
    #print("\n SX1: ",SX1)
    #print("\n SX2: ",SX2)
    plot(iterl, FArray, linetype = :path, title = "Line Plot", xlabel = "n_iter", ylabel = "objective function", color = :red)
    plot(SX1,SX2, linetype = :path, title = "Line Plot", xlabel = "sx1", ylabel = "sx2", color = :red)
    plot(lsx,lxz, linetype = :path, title = "Line Plot", xlabel = "x1", ylabel = "x2", color = :red)