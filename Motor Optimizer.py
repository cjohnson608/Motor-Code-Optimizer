import numpy as np
from scipy.optimize import minimize
numGrains = 5.0
odia = 6.0
dia = 5.0
height = 4.0
A_t = 1.5
x0 = np.array([numGrains, odia, dia, height, A_t])

def pseudo_prop_sim(x0):
    targetName = 'Impulse'
    number, odia, dia, height, A_t = x0
    number = round(number)

    m_flux,peak_pressure,w_grain,peak_thrust, tot_impulse,number,odia,dia,height,A_t,thrustVec = prop_sim(number,odia,dia,height,A_t)

    if targetName == 'Thrust':
        targetVal = peak_thrust
    elif targetName == 'Impulse':
        targetVal = tot_impulse
        
    if (type(x0) != float) or (type(targetVal) != float):
        targetVal = 0
        print("Internal error: imaginary values. Try a new start point.")
    return targetVal * -1


def prop_sim(number, odia, dia, height, A_t):
    a = 0.027 #coeff of pressure
    n = 0.3 #burning rate exponent (dimensionless)
    gamma = 1.2 #specific heat ratio (dimensionless)
    M = 2 #mach
    R = 2000 #gas constant
    T = 4700 #temp R
    c_star = 4890 #ft/sec
    rho_prop = 0.00176 #slug/in^3
    p_a = 14.7 #psi
    num = number
    inhib = 0

    w_grain = 0.5*(odia-dia)
    dia_vec = dia*np.ones(num)
    h_vec = height*np.ones(num)

    #isentropic relationships
    exp_r = 1.884
    A_exit = exp_r*A_t #exit area from expansion ratio
    t_e = T*((1+((gamma-1)/2)*M**2)**-1) #adjust temp

    #initialize
    Rad = 0.5*odia
    Pressure = np.array([])
    Thrust = np.array([])
    time = np.array([])
    Burn_Rate = np.array([])
    Mass_flow = np.array([])

    #find distance steps

    odia_vec = odia*np.ones(num)
    diff_vec = odia_vec - dia_vec
    max_dis = np.max(diff_vec)
    d_step = max_dis/2000
    d_step_vec = d_step*np.ones(1000)

    #mass flux calc

    A_port_start = (np.pi*(5*dia)**2)
    A_port = np.array([A_port_start])
    port_dia = np.array([dia])

    test = np.append(port_dia, 200)

    for q in range(0,d_step_vec.size-1):
        port_dia = np.append(port_dia, port_dia[q] + 2*(d_step))
        A_port = np.append(A_port, (np.pi*(5*port_dia[q])**2))

    Total_A_b_vec = np.array([])
    m = 0 #not needed !!!!
    for y in range(0,d_step_vec.size):
        grains = np.array([])
        for w in range(0, num-1):
            dia = dia_vec[w]
            h = h_vec[w]
            rad = 0.5*dia
            a_inside = np.pi*dia*h
            a_circle = (np.pi*Rad**2)-(np.pi*rad**2)
            a_circle_total = a_circle*(2-inhib)
            A_grain = a_inside+a_circle_total
            diff = odia-dia
            if diff<=0:
                A_grain = 0;
            grains = np.append(grains, A_grain)
        Burn_Area = sum(grains)
        Total_A_b_vec = np.append(Total_A_b_vec, Burn_Area)
        dia_vec = dia_vec + 2*d_step
        if inhib ==0:
            h_vec = h_vec - 2*d_step
        m = m+1
    
    i = 0 # Not needed!!!

    for c in Total_A_b_vec:
        A_b = c
        K_n = A_b/A_t
        P_c = (K_n*c_star*a*rho_prop)**(1/(1-n))
        #Thrust
        #Calc m_dpt
        r_b = a*P_c**n
        m_dot = rho_prop*A_b*r_b
        #Calc exit pressure
        ratio_p = (1+((gamma-1)/2)*(M**2))**(-gamma/(gamma-1))
        p_exit = P_c*ratio_p
        #calc exit velocity
        v_e = M*(np.sqrt(gamma*R*t_e))
        F = (m_dot*v_e)+A_exit*(p_exit-p_a)
        Pressure = np.append(Pressure, P_c)
        Thrust = np.append(Thrust, F)
        Mass_flow = np.append(Mass_flow, m_dot)

        time_step = d_step/r_b
        if time.size == 0:
            time = np.append(time, time_step)
        else:
            time = np.append(time, time[time.size-1] + time_step)
        
        i = i +1

    m_flux_vec = Mass_flow/A_port

    m_flux = np.max(m_flux_vec)*32.17
    peak_thrust = np.max(Thrust)
    peak_pressure = max(Pressure)
    if time.size == 0:
        tot_impulse = 0
    else:
        tot_impulse = np.trapz(Thrust, time) #CHEK THIS LINE, REVERSED
    #cell array stuff here?
    thrustVec = Thrust
    print("total impulse: ")
    print(tot_impulse)
    return (m_flux,peak_pressure,w_grain,peak_thrust, tot_impulse,number,odia,dia,height,A_t,thrustVec)

bnds = ((1, 15), (odia-1, odia+1), (0.01, odia), (3, np.Infinity), (0.5, 2))

m_flux,peak_pressure,w_grain,peak_thrust, tot_impulse,number,odia,dia,height,A_t,thrustVec = prop_sim(int(x0[0]),float(x0[1]),float(x0[2]),float(x0[3]),float(x0[4]))
def c1(x0):
    return x0[1]-x0[2]
def c2(x0):
    return 72 - x0[0]*x0[3]
def c3(peak_pressure):
    return peak_pressure - 300
def c4(peak_pressure):
    return 500 - peak_pressure
def c5(m_flux):
    return 1.2 - m_flux
def c6(x0):
    return 3.141592653589793*((x0[2]/2)**2)-x0[4]**2
def c7(peak_thrust):
    return peak_thrust

cons = [{'type':'ineq', 'fun': c1},
        {'type':'ineq', 'fun': c2},
        {'type':'ineq', 'fun': c3},
        {'type':'ineq', 'fun': c4},
        {'type':'ineq', 'fun': c5},
        {'type':'eq', 'fun': c6},
        {'type':'ineq', 'fun': c7}
        ]
res = minimize(pseudo_prop_sim, x0, args=(), method='SLSQP', jac=None, bounds=bnds, constraints=cons, tol=None, callback=None, options={'maxiter': 100, 'ftol': 1e-06, 'iprint': 1, 'disp': True, 'eps': 1.4901161193847656e-08, 'finite_diff_rel_step': None})
#res = minimize(pseudo_prop_sim, x0, method='SLSQP', jac = None, bounds = bnds, constraints = cons, options={'xtol': 1e-3, 'disp': True})
#res = minimize(pseudo_prop_sim, x0, args=(), method='COBYLA', constraints=cons, tol=None, callback=None, options={'rhobeg': 1.0, 'maxiter': 1000, 'disp': True, 'catol': 0.0002})
print("Solution array: ")
print(res.x)
print(prop_sim(5, 6, 2.4775, 9.6375, 2))