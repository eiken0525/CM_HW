import math
import matplotlib.pyplot as plt

def plot_array(x_axis, y_axis, title):
    plt.figure(figsize=(10, 6))
    plt.plot(x_axis, y_axis, 'o', label='Approximated Solution')
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(f"P{title[-3:]}.png", transparent=True)

def adams_variable_step_size_predictor_corrector(f, t, y, b, h_min, h_max, tol, N):
    def rk4(f, t, y, h):
        k1 = h*f(t, y)
        k2 = h*f(t + h/2, y + k1/2)
        k3 = h*f(t + h/2, y + k2/2)
        k4 = h*f(t + h, y + k3)
        t=t+h
        return t,y + (k1 + 2*k2 + 2*k3 + k4)/6
    
    AB_correctors = [
        [2, 3, -1],
        [12, 23, -16, 5],
        [24, 55, -59, 37, -9],
        [720, 1901, -2774, 2616, -1274, 251]
    ]
    AM_correctors=[
        [2, 1, 1],
        [12, 5, 8, -1],
        [24, 9, 19, -5,1],
        [720, 251, 646, -264, 106,-19]   
    ]

    corrector = AB_correctors[N-2]
    time = [t]
    y_values = [y]
    result=[]
    last=False
    h=h_max

    for _i in range(1, N):
        t_i,y_i=rk4(f, time[_i-1], y_values[_i-1],h)
        y_values.append(y_i)
        time.append(t_i)
    nflag=True
    i=N
    t=time[len(time)-1]+h

    while h != 0:
        fix = 0
    
        for j in range(N):
            fix += corrector[j+1] * f(time[i-j-1], y_values[i-j-1])
        wp = y_values[i-1] + h/corrector[0]*fix

        
        coeffs = AM_correctors[N-2]
        fix = coeffs[1]*f(t, wp)
    
        for j in range(len(coeffs)-2):
            fix += coeffs[j+2] * f(time[i-j-1], y_values[i-j-1])
        wc = y_values[i-1] + h/coeffs[0]*fix

        sigma=19*abs(wc-wp)/(270*h)

        if(sigma<=tol):
            y_values.append(wc)
            time.append(t)
            if nflag:
                for j in range(N):
                    result.append([time[i-N+j],y_values[i-N+j]])
            else:
                result.append([time[i],y_values[i]])
            if last:
                # print('last')
                result.append([time[i],y_values[i]])
                break
            else:
                i=i+1

                nflag=False
                if sigma<=(0.1*tol) or time[i-1]+h>b:
                    q=(tol/(2*sigma)) ** (0.25)
                    if q>4:
                        h=4*h
                    else: 
                        h=q*h
                        
                    if h>h_max:
                        h=h_max
                        
                    if time[i-1]+4*h>=b:
                        h=(b-time[i-1])/4
                        last=True
                    
                    for _i in range(-1, N-2):
                        t_i,y_i=rk4(f, time[i+_i], y_values[i+_i],h)
                        y_values.append(y_i)
                        time.append(t_i)
                    nflag=True
                    
                elif time[i-1]+4*h>=b:
                    h=(b-time[i-1])/4
                    last=True

                    for _i in range(-1, N-2):
                        t_i,y_i=rk4(f, time[i+_i], y_values[i+_i],h)
                        y_values.append(y_i)
                        time.append(t_i)
                    nflag=True
                    i=i+N-1
        else:
            q=(tol/(2*sigma))**0.25

            if(q<0.1):
                h=0.1*h
            else: 
                h=q*h
            
            if h<h_min:
                break
            else :
                if nflag:
                    i=i-(N-1)

                for _i in range(-1, N-2):
                        t_i,y_i=rk4(f, time[i+_i], y_values[i+_i],h)
                        if((i+_i+1)>=len(y_values)):
                            y_values.append(y_i)
                            time.append(t_i)
                            
                        else:
                            y_values[i+_i+1]=y_i
                            time[i+_i+1]=t_i
                nflag=True
                i=i+N-1
        t=time[len(time)-1]+h
    t_values, approx_soln_list = [result[i][0] for i in range(len(result))], [result[i][1] for i in range(len(result))]
    return t_values, approx_soln_list


def f_a(t, y):
    return math.sin(t) + math.exp(-t)

def f_b(t, y):
    return - t * y + 4 * t / y

t_values_a, approximation_a = adams_variable_step_size_predictor_corrector(f_a, 0, 0, 1, 0.01, 0.2, 1e-4, 4)

plot_array(t_values_a, approximation_a, "Approximation by the Adams Variable Step-Size Predictor-Corrector Algorithm: 13a")

t_values_b, approximation_b = adams_variable_step_size_predictor_corrector(f_b, 0, 1, 1, 0.01, 0.2, 1e-4, 4)

plot_array(t_values_b, approximation_b, "Approximation by the Adams Variable Step-Size Predictor-Corrector Algorithm: 13b")