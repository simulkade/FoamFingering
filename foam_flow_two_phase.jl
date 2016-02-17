"""
foam flow in porous media: scaling the relperm of gas phase
Written by AA Eftekhari
Delft, November 2015

Usage:
```
using JFVM, ProgressMeter
m = 80 # number of cells in x direction
n = 50 # number of cells in y direction
W = 10 # width
H = 1 # height
x=[linspace(0,0.2,100); linspace(0.21, 1.0, 10)]
y=collect(linspace(0, H, n))
meshvar = createMesh1D(x)
foam_stars(meshvar)
```
"""
function foam_stars(meshvar; muw=0.001, mug=2e-5, perm_ave=1e-12,
  poros_ave=0.2, fmmob=25000.0, epdry=10000.0, fmdry=0.29,
  swc=0.1, sgr=0.05, krg0=0.9, ng=1.8, krw0=0.2, nw=4.2, sw_init=1.0,
  t_final = 6.145, dt0 = 0.01)
# Physical properties and rel-perms
sws(sw::Real)=((sw>swc)*(sw<1-sgr)*(sw-swc)/(1-sgr-swc)+(sw>=1-sgr)*1.0)
sws(sw::Array{Float64})=((sw.>swc).*(sw.<1-sgr).*(sw-swc)/(1-sgr-swc)+(sw.>=1-sgr).*ones(size(sw)))
kr(sw)=(krg0*(1-sws(sw)).^ng)
fm(sw)=(1+fmmob*(0.5+atan(epdry.*(sw-fmdry))/π))
krg(sw)=(kr(sw)./fm(sw))
krw(sw)=(krw0*sws(sw).^nw)
dkrwdsw(sw)=(nw*krw0*(1/(1-sgr-swc))*sws(sw).^(nw-1))
dkrdsw(sw)=((krg0*ng*(1-sws(sw)).^(ng-1))/(-swc-sgr+1))
dfmdsw(sw)=(((epdry*fmmob)./(π*(epdry^2*(sw-fmdry).^2+1))))
dkrgdsw(sw)=((dkrdsw(sw).*fm(sw)-dfmdsw(sw).*kr(sw))./fm(sw).^2)
fw(sw)=((krw(sw)/muw)./(krw(sw)/muw+krg(sw)/mug))
dfw(sw)=((dkrwdsw(sw)/muw.*(krw(sw)/muw+krg(sw)/mug)-
    (dkrwdsw(sw)/muw+dkrgdsw(sw)/mug).*krw(sw)/muw)./
    (krg(sw)/mug+krw(sw)/muw).^2)
# Assign the properties to the cells
poros = createCellVariable(meshvar, poros_ave)
perm = createCellVariable(meshvar, perm_ave)
mu_gas = createCellVariable(meshvar, mug) # gas viscosity
mu_water = createCellVariable(meshvar, muw) # water viscosity
Lg_ave = harmonicMean(perm./mu_gas)
Lw_ave = harmonicMean(perm./mu_water)
s = createCellVariable(meshvar, 0) # gas source term
# BC
BCp = createBC(meshvar) # all Neumann BC for pressure
BCs = createBC(meshvar) # saturation BC
if meshvar.dimension<2
  BCp.left.a[:] = krg(0)*Lg_ave.xvalue[1]
elseif 2<=meshvar.dimension<3
  BCp.left.a[:] = krg(0)*Lg_ave.xvalue[1,:]
  BCp.top.periodic=true
  BCp.bottom.periodic=true
  BCs.top.periodic=true
  BCs.bottom.periodic=true
end
BCp.left.b[:] = 0
u = 1e-3
BCp.left.c[:] = -u
BCp.right.a[:]=0
BCp.right.b[:] = 1
BCp.right.c[:] = 1e5
BCs.left.a[:] = 0
BCs.left.b[:] = 1
BCs.left.c[:] = 0.0
# initial condition
sw0 = createCellVariable(meshvar, sw_init, BCs) # initial gas saturation
s = createCellVariable(meshvar, 0) # gas source term
p0 = createCellVariable(meshvar, 1e5, BCp) # [Pa] initial pressure
sw_old = copyCell(sw0)
p_old = copyCell(p0)
sw = copyCell(sw_old)
p = copyCell(p_old)
pgrad = gradientTerm(p)
# solver setting
dt=dt0 # [s] time step
t=0
eps_sw = 1e-5
eps_p = 1e-2
# Explicit M and RHS
RHSs = constantSourceTerm(s) # explicit source term to be added to the rhs
(BCMp, BCRHSp) = boundaryConditionTerm(BCp)
pr=Progress(100, 0.1, "Time loop:", 50)
while t<t_final
    error1 = 1e5
    error2 = 1e5
    loop_count =0
    while true
        loop_count=loop_count+1
        if loop_count>10
            p=copyCell(p_old)
            sw=copyCell(sw_old)
            dt=dt/5
            break
        end
        if (error1<=eps_sw)&&(error2<=eps_p)
            t=t+dt
            dt=dt0
            sw_old = copyCell(sw)
            p_old = copyCell(p)
            #println(t)
            #plot(log(L_ave.xvalue')) #shading interp drawnow
            #figure(1)  pcolor(1-sw.value(2:m+1,2:n+1)') shading interp colorbar drawnow
            #figure(1) visualizeCells(meshvar, sw.value) drawnow
#             figure(1)
#             subplot(2,1,1)plot(sw_ave.xvalue, 'o')
#             ylabel('S_w')
#             subplot(2,1,2) semilogy(L_ave.xvalue)
#             xlabel('x')
#             ylabel('total mobility')
#             drawnow
#             plot(t, p.value(2), 'o') drawnow
#             plot(L_ave.xvalue) drawnow
            break
        end
        # step 1) calculate the average values
#         sg_ave = arithmeticMean(meshvar, sg.value)
        sw_ave = upwindMean(sw, -pgrad)
        Lg = Lg_ave.*faceEval(krg,sw_ave)
        Lw = Lw_ave.*faceEval(krw,sw_ave)

        # step 2) calculate the pressure profile
        L_ave = Lw+Lg
        Meq = diffusionTerm(L_ave)


        # solve the linear system of equations and reshape the result
        Mp = Meq + BCMp
        RHSp = BCRHSp - RHSs # the whole continuity is multiplied by a minus sign
        p_new = solvePDE(meshvar, Mp, RHSp)
#         P = Mp\RHSp
#         p.value = reshape(full(P), m+2, n+2)

        pgrad = gradientTerm(p_new)
        error1=1e5
#         while error1>eps_sw
        for i = 1:3
            sw_ave = upwindMean(sw, -pgrad)
    #         Lg = Lg_ave.*krg(sw_ave)
            Lw = Lw_ave.*faceEval(krw,sw_ave)
    #         sw_ave = upwindMean(meshvar, -pgrad, sw.value)
        #     sw_ave = arithmeticMean(meshvar, sw.value)
    #         Lw = Lw_ave.*krw(sw_ave)
        #     Lo = Lo_ave.*kro(sw_ave)
            # step 3) calculate the new value of sw
            (Mtrans, RHStrans) = transientTerm(sw_old, dt, poros)
            u = -faceEval(dkrwdsw,sw_ave).*Lw_ave.*pgrad
            Mconv = convectionUpwindTerm(u)
    #       Mconv = convectionTerm(meshvar, u)
            facevar = (-Lw+faceEval(dkrwdsw,sw_ave).*Lw_ave.*sw_ave).*pgrad
            RHSdiv = divergenceTerm(facevar)
            (BCM, BCRHS) = boundaryConditionTerm(BCs)

            # construct the linear system
            M = Mtrans+Mconv+BCM
            RHS = RHStrans+BCRHS-RHSdiv+RHSs
            sw_new = solvePDE(meshvar, M, RHS)
    #         SG = M\RHS
    #         sg_new = reshape(full(SG), m+2, n+2)
            error1 = sum(abs(sw_new.value[:]-sw.value[:]))
#             sw.value = 0.3*sw_new+0.7*sw.value
            sw=copyCell(sw_new)
        end
        error2 = maximum(abs(p_new.value-p.value)./p_new.value)
#         p.value = 0.1*p.value+0.9*p_new
        p = copyCell(p_new)
    end
    update!(pr, round(Int,t*100/t_final+1.0)) # progress bar
end
return sw
end #function
