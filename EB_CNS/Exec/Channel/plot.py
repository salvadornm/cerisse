# Monitor the total quantities extracted from the screen output
import matplotlib.pyplot as plt
import re
import numpy as np
import yt
import os
import argparse

plt.style.use("yt.default")
plt.rcParams["font.size"] = 12

try:
    default_fname = sorted([f for f in os.listdir(".") if f.startswith("plt")])[-1]
except:
    default_fname = None

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fname", default=default_fname, help="pltfile to read")
parser.add_argument("-s", "--scn", default=None, help="screen file to read")
parser.add_argument("-r", "--Reb", default=1.7e4, type=float, help="bulk Reynolds")
parser.add_argument("-m", "--Mb", default=1.5, type=float, help="bulk Mach")
parser.add_argument("-w", "--wm", default=False, type=bool, help="applied wall model?")
args = parser.parse_args()
print(args)

# Problem settings
Tw = 500.0
H = 0.6845
muw = 2.67e-4
ub = 44822.0654 * args.Mb
rhob = muw/ub/H * args.Reb

# 1. Read screen output
if args.scn is not None:
    scnname = args.scn
    print("Reading screen output from " + scnname + " ...")
    quantities = {}
    # Define a regex pattern to match the lines with quantities
    pattern = re.compile(r"Total (\w+[-]?\w+)\s+=\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
    with open(scnname, 'r') as file: # Open the file and read it line by line
        for line in file:
            # Use the regex pattern to search for a match in the line
            match = pattern.search(line)
            if match:
                # Extract the quantity name and value
                quantity_name = match.group(1)
                quantity_value = float(match.group(2))
                # Store them in the dictionary
                if quantity_name in quantities:
                    quantities[quantity_name].append(quantity_value)
                else:
                    quantities[quantity_name] = [quantity_value]

    plt.figure(figsize=(4, 8), constrained_layout=True)
    for i, (key, value) in enumerate(quantities.items()):
        plt.subplot(len(quantities.keys()),1,i+1)
        plt.plot(value)
        plt.ylabel(key)
        plt.xlim(0, len(value))
    plt.xlabel('step')
    print(" >>> last step: " + str(len(value)))
    plt.savefig("plot1.png")
    print(" >>> Screen output saved to plot1.png")

if args.fname is not None:
    fname = args.fname
    print("Opening pltfile " + fname + " ...")

    # 2. Plot statistics
    ds = yt.load(fname)
    ad = ds.all_data()
    y = ad[("y")].to_ndarray()
    rho = ad[("density")].to_ndarray()
    u = ad[("x_velocity")].to_ndarray()
    v = ad[("y_velocity")].to_ndarray()
    w = ad[("z_velocity")].to_ndarray()
    T = ad[("Temp")].to_ndarray()
    mu = ad[("viscosity")].to_ndarray()
    lam = ad[("conductivity")].to_ndarray()
    cp = np.mean(ad[("cp")].to_ndarray())

    # averaged conditionally on y
    y_unique = np.unique(y) # find unique y values
    mean = {key: np.zeros_like(y_unique) for key in ["rho", "u", "v", "w", "T", "mu", "lambda"]}
    var = {key: np.zeros_like(y_unique) for key in ["u", "v", "w", "uv", "T"]}
    for i, y_val in enumerate(y_unique):
      mean["rho"][i] = np.mean(rho[y == y_val])
      mean["u"][i] = np.mean(u[y == y_val])
      mean["v"][i] = np.mean(v[y == y_val])
      mean["w"][i] = np.mean(w[y == y_val])
      mean["T"][i] = np.mean(T[y == y_val])
      mean["mu"][i] = np.mean(mu[y == y_val])
      mean["lambda"][i] = np.mean(lam[y == y_val])
      var["u"][i] = np.mean(u[y == y_val]**2)
      var["v"][i] = np.mean(v[y == y_val]**2)
      var["w"][i] = np.mean(w[y == y_val]**2)
      var["uv"][i] += np.mean(u[y == y_val]*v[y == y_val])
      var["T"][i] = np.mean(T[y == y_val]**2)
    var["u"] -= mean["u"]**2
    var["v"] -= mean["v"]**2
    var["w"] -= mean["w"]**2
    var["uv"] -= mean["u"]*mean["v"]
    var["T"] -= mean["T"]**2
    y = y_unique / H
    print(f' >>> max(u/ub) = {np.max(mean["u"])/ub}; max(T/Tw) = {np.max(mean["T"])/Tw}; max(rho/rhob) = {np.max(mean["rho"])/rhob}')

    ref_data = np.loadtxt("./coleman1995.dat", skiprows=2, delimiter=',')
    linestyle = '.-'
    plt.figure(figsize=(5, 6), constrained_layout=True)
    plt.subplot(3,1,1)
    plt.plot(y, mean["u"]/ub, linestyle, label=r"$\langle u\rangle/u_b$", color="C0")
    plt.plot(y, mean["v"]/ub, linestyle, label=r"$\langle v\rangle/u_b$", color="C1")
    plt.plot(y, mean["w"]/ub, linestyle, label=r"$\langle w\rangle/u_b$", color="C2")
    plt.plot(y, mean["T"]/Tw, linestyle, label=r"$\langle T\rangle/T_w$", color="C3")
    plt.plot(y, mean["rho"]/rhob, linestyle, label=r"$\langle\rho\rangle/\rho_b$", color="C4")
    plt.plot(ref_data[:,0], ref_data[:,4], '--', color="C0") # u
    plt.plot(ref_data[:,0], ref_data[:,3], '--', color="C3") # T
    plt.plot(ref_data[:,0], ref_data[:,1], '--', color="C4") # rho
    plt.xlim([-1, 1])
    plt.ylim([0, None])
    plt.xticks([])
    plt.legend()

    plt.subplot(3,1,2)
    plt.plot(y, var["u"]/ub**2, linestyle, label=r"$\langle u'u'\rangle/u_b^2$", color="C0")
    plt.plot(y, var["v"]/ub**2, linestyle, label=r"$\langle v'v'\rangle/u_b^2$", color="C1")
    plt.plot(y, var["w"]/ub**2, linestyle, label=r"$\langle w'w'\rangle/u_b^2$", color="C2")
    #plt.plot(y, var["uv"]/ub**2, linestyle, label=r"$\langle u'v'\rangle/u_b^2$", color="C5")
    #plt.plot(y, var["T"]/Tw**2, linestyle, label=r"$\langle T'T'\rangle/T_w^2$", color="C3")
    plt.plot(ref_data[:,0], ref_data[:,16], '--', color="C0") # u'u'
    plt.plot(ref_data[:,0], ref_data[:,17], '--', color="C1") # v'v'
    plt.plot(ref_data[:,0], ref_data[:,18], '--', color="C2") # w'w'
    #plt.plot(ref_data[:,0], ref_data[:,19], '--', color="C5") # u'v'
    #plt.plot(ref_data[:,0], ref_data[:,20], '--', color="C3") # T'T'
    plt.xlim([-1, 1])
    plt.ylim([0, None])
    plt.xticks([])
    plt.legend()

    plt.subplot(3,1,3)
    plt.plot(y, var["uv"]/ub**2, linestyle, label=r"$\langle u'v'\rangle/u_b^2$", color="C5")
    plt.plot(y, var["T"]/Tw**2, linestyle, label=r"$\langle T'T'\rangle/T_w^2$", color="C3")
    plt.plot(ref_data[:,0], ref_data[:,19], '--', color="C5") # u'v'
    plt.plot(ref_data[:,0], ref_data[:,20], '--', color="C3") # T'T'
    plt.xlim([-1, 1])
    plt.legend()
    plt.savefig("plot2.png")
    print(" >>> Statistics saved to plot2.png")

    # 3. Plot U+ vs y+
    if not args.wm:
        tauw = mean["mu"][0] * (mean["u"][0] - 0) / ((y[0]+1)*H)
        qw = mean["lambda"][0] * (mean["T"][0] - Tw) / ((y[0]+1)*H)
        rhow = mean["rho"][0]
        utau = np.sqrt(tauw / rhow)
        yplus = (-np.abs(y) + 1) * H * utau / (mean["mu"][0] / rhow)
    else:        
        def viscosity(T):
          return muw*(T/Tw)**0.7

        def conductivity(T):
          mu = viscosity(T)
          Pr = 0.72
          return cp / Pr * mu
  
        def eqwm(hwm, rhowm, uwm, Twm, Tw, couple_T=True, N=100):
            # Grid param
            stretch = 1.1 # stretch factor
            # Modelling constants
            kappa = 0.41
            Aplus = 17.
            Prt = 0.9
            # Wall condition
            uw = 0.0

            y = hwm*(stretch**np.arange(N+1)-1)/((stretch**N+stretch**(N-1))/2-1)
            y = (y[1:]+y[:-1])/2
            u = y / hwm * (uwm-uw) + uw # initial u
            T = y / hwm * (Twm-Tw) + Tw # initial T

            res_u = 1e10
            res_T = 0
            while res_u > 0.001 * uwm or res_T > 0.001 * Twm: # 0.1% error
              mu = viscosity(T)
              tauw = mu[0]*(u[0]-uw)/y[0]
              rho = rhowm * Twm / T
              yplus = y * np.sqrt(tauw / rho[0]) / (mu[0] / rho[0])
              mut = kappa * y * np.sqrt(rho * tauw) * (1. - np.exp(-yplus / Aplus))**2

              # Solve for u
              A = np.zeros((N, N))
              b = np.zeros((N,))
              # Wall BC
              mu_lo = mu[0]+mut[0]
              mu_hi = 0.5*(mu[0]+mut[0]+mu[1]+mut[1])
              A[0,0] = mu_lo*(-1./y[0]) - mu_hi*(1./(y[1]-y[0]))
              A[0,1] = mu_hi*(1./(y[1]-y[0]))
              b[0] = -mu_lo*uw/y[0]
              # Freestream BC
              A[-1,-1] = 1.
              b[-1] = uwm
              # In the middle
              for i in range(1,N-1):
                mu_lo = 0.5*(mu[i]+mut[i]+mu[i-1]+mut[i-1])
                mu_hi = 0.5*(mu[i]+mut[i]+mu[i+1]+mut[i+1])
                A[i,i-1] = mu_lo/(y[i]-y[i-1])
                A[i,i] = -mu_lo/(y[i]-y[i-1]) - mu_hi/(y[i+1]-y[i])
                A[i,i+1] = mu_hi/(y[i+1]-y[i])
              u2 = np.linalg.solve(A, b)
              res_u = np.sqrt(np.mean((u2-u)**2))
              u = u2

              if couple_T:
                # Solve for T
                lam = conductivity(T)
                C = np.zeros((N, N))
                d = np.zeros((N,))
                if Tw > 0:
                  # Wall BC (isothermal)
                  mu_lo = mu[0]+mut[0]
                  mu_hi = 0.5*(mu[0]+mut[0]+mu[1]+mut[1])
                  lam_lo = (lam[0]+cp*mut[0]/Prt)
                  lam_hi = 0.5*(lam[0]+cp*mut[0]/Prt+lam[1]+cp*mut[1]/Prt)
                  C[0,0] = lam_lo*(-1./y[0]) - lam_hi*(1./(y[1]-y[0]))
                  C[0,1] = lam_hi*(1./(y[1]-y[0]))
                  d[0] = -lam_lo*Tw/y[0] + mu_lo*(uw*u[0]/y[0]) - mu_hi*((u[0]+u[1])/2*(u[1]-u[0])/(y[1]-y[0]))
                else:
                  # Wall BC (adiabatic)
                  mu_lo = mu[0]+mut[0]
                  mu_hi = 0.5*(mu[0]+mut[0]+mu[1]+mut[1])
                  lam_lo = (lam[0]+cp*mut[0]/Prt)
                  lam_hi = 0.5*(lam[0]+cp*mut[0]/Prt+lam[1]+cp*mut[1]/Prt)
                  C[0,0] = - lam_hi*(1./(y[1]-y[0]))
                  C[0,1] = lam_hi*(1./(y[1]-y[0]))
                  d[0] = mu_lo*(uw*u[0]/y[0]) - mu_hi*((u[0]+u[1])/2*(u[1]-u[0])/(y[1]-y[0]))
                # Freestream BC
                C[-1,-1] = 1.
                d[-1] = Twm
                # In the middle
                for i in range(1,N-1):
                  mu_lo = 0.5*(mu[i]+mut[i]+mu[i-1]+mut[i-1])
                  mu_hi = 0.5*(mu[i]+mut[i]+mu[i+1]+mut[i+1])
                  lam_lo = 0.5*(lam[i-1]+cp*mut[i-1]/Prt+lam[i]+cp*mut[i]/Prt)
                  lam_hi = 0.5*(lam[i]+cp*mut[i]/Prt+lam[i+1]+cp*mut[i+1]/Prt)
                  C[i,i-1] = lam_lo/(y[i]-y[i-1])
                  C[i,i] = -lam_lo/(y[i]-y[i-1]) - lam_hi/(y[i+1]-y[i])
                  C[i,i+1] = lam_hi/(y[i+1]-y[i])
                  d[i] = mu_lo*((u[i-1]+u[i])/2*(u[i]-u[i-1])/(y[i]-y[i-1])) - mu_hi*((u[i]+u[i+1])/2*(u[i+1]-u[i])/(y[i+1]-y[i]))
                T2 = np.linalg.solve(C, d)
                res_T = np.sqrt(np.mean((T2-T)**2))
                T = T2

            mu0 = viscosity(T[0])
            tauw = mu0*(u[0]-uw)/y[0]
            lam0 = conductivity(T[0])
            qw = lam0*(T[0]-Tw)/y[0]
            return tauw, qw, y, u, T
        
        wmpt = 1
        tauw, qw, ywm, uwm, Twm = eqwm((y[wmpt]+1)*H, mean["rho"][wmpt], mean["u"][wmpt], mean["T"][wmpt], Tw)
        rhow = mean["rho"][wmpt] * mean["T"][wmpt] / Tw
        utau = np.sqrt(tauw / rhow)
        yplus = (-np.abs(y) + 1) * H * utau / (muw / rhow)
    print(" >>> Wall model =", args.wm, ": utau =", utau, ", -Bq =", qw/(mean["rho"][0]*cp*utau*Tw), ", tauw =", tauw)
    uplus = mean["u"] / utau

    uplusvd = np.zeros_like(yplus)
    uplusvd[0] = uplus[0]
    uplusvd[-1] = uplus[-1]
    for i in range(1, len(uplus)//2+1):
      uplusvd[i] = ((mean["rho"][i]/rhow)**0.5 * (uplus[i]-uplus[i-1])) + uplusvd[i-1]
      uplusvd[-i] = ((mean["rho"][-i]/rhow)**0.5 * (uplus[-i]-uplus[-i+1])) + uplusvd[-i+1]

    def myplot(*args, **kwargs):
      # plt.plot(*args, **kwargs)
      plt.semilogx(*args, **kwargs)

    dnsdata = np.loadtxt("./modesti2016.dat", skiprows=2, delimiter=',')

    plt.figure(layout="constrained")
    plt.subplot(1,2,1)
    myplot(dnsdata[:,0], dnsdata[:,1:], 'o', color="grey", label=["DNS",None,None,None,None])
    myplot(yplus, uplusvd, '.-', label="$U^+_{VD}$")
    myplot(yplus, var["u"]/utau**2*mean["rho"]/rhow, '.-', label=r"$\tau_{VD, 11}/u_{\tau}^2$")
    myplot(yplus, var["v"]/utau**2*mean["rho"]/rhow, '.-', label=r"$\tau_{VD, 22}/u_{\tau}^2$")
    myplot(yplus, var["w"]/utau**2*mean["rho"]/rhow, '.-', label=r"$\tau_{VD, 33}/u_{\tau}^2$")
    myplot(yplus, -np.abs(var["uv"])/utau**2*mean["rho"]/rhow, '.-', label=r"$\tau_{VD, 12}/u_{\tau}^2$")
    ypp = np.linspace(1, 10.8, 20)
    myplot(ypp, ypp, 'k--')
    ypp = np.linspace(11, yplus[len(yplus)//2], 20)
    myplot(ypp, np.log(ypp) / 0.41 + 5, 'k--')

    if args.wm:
        ypluswm = ywm * utau / (muw / rhow)
        upluswm = uwm / utau
        uplusvdwm = np.zeros_like(ypluswm)
        uplusvdwm[0] = ypluswm[0] #upluswm[0]
        for i in range(1, len(upluswm)):
            rhoi = mean["rho"][wmpt] * mean["T"][wmpt] / Twm[i]
            uplusvdwm[i] = ((rhoi/rhow)**0.5 * (upluswm[i]-upluswm[i-1])) + uplusvdwm[i-1]
        myplot(ypluswm, upluswm, '--', color="C0")

    plt.xlim([1, 1000])
    plt.xlabel("$y^+$")
    #plt.ylabel("$U^+_{VD}$")
    plt.legend()

    plt.subplot(1,2,2)
    ue = mean["u"][len(mean["u"])//2]
    Te = mean["T"][len(mean["T"])//2]
    plt.plot(mean["u"]/ue, mean["T"]/Tw, '.-')
    rg = 2*cp*(Tw-Te)/ue**2 - 2*0.72*(-qw)/ue/tauw # take Pr=0.72, and qw as -qw
    Trg = Te + rg*ue**2/2/cp
    plt.plot(mean["u"]/ue, 1 + (Trg-Tw)/Tw*(mean["u"]/ue) + (Te-Trg)/Tw*(mean["u"]/ue)**2, 'k--', label="Zhang")
    plt.xlabel("$u/u_e$")
    plt.ylabel("$T/T_w$")
    plt.legend()
    plt.savefig("plot3.png")
    print(" >>> Non-dimensional profiles saved to plot3.png")

    # 4. Slice plot
    slc = yt.SlicePlot(ds, "z", "x_velocity", center=(4.107, 0.0, 0.0))
    slc.set_log("x_velocity", False)
    slc.save("plot4.png")
    slc = yt.SlicePlot(ds, "z", "Temp", center=(4.107, 0.0, 0.0))
    slc.set_log("Temp", False)
    slc.save("plot5.png")

    # (SF only) 5. ksgs
    try:
        y = ad[("y")].to_ndarray()
        y_unique = np.unique(y)
        R11_raw = ad[("R11")].to_ndarray()
        R22_raw = ad[("R22")].to_ndarray()
        R33_raw = ad[("R33")].to_ndarray()
        R11 = np.zeros_like(yplus)
        R22 = np.zeros_like(yplus)
        R33 = np.zeros_like(yplus)
        for i, y_val in enumerate(y_unique):
            R11[i] = np.mean(R11_raw[y == y_val])
            R22[i] = np.mean(R22_raw[y == y_val])
            R33[i] = np.mean(R33_raw[y == y_val])
        plt.figure(layout="constrained")
        myplot(yplus, (R11+R22+R33)/2, '.-', label=r"$k^{sgs}$")
        myplot(yplus, R11, '.-', label=r"$\tau_{11}^{sgs}$")
        myplot(yplus, R22, '.-', label=r"$\tau_{22}^{sgs}$")
        myplot(yplus, R33, '.-', label=r"$\tau_{33}^{sgs}$")
        myplot(yplus, (var["u"]+var["v"]+var["w"])/2/10, '.-', label=r"$k_t/10$")
        plt.legend()
        plt.savefig("plot6.png")
        print(" >>> ksgs profiles saved to plot6.png")
    except:
        ksgs = 0

print("DONE")
