import sys
import matplotlib.pyplot as plt
import numpy as np

print("Solves dx(t)/dt = -death * x(t) + birth * x(t - lag) / (1 + x(t - lag) / qm)")
print("Using:")
print(" - Explicit Euler: x(t + dt) = x(t) + dt * (-death * x(t) + birth * x(t - lag) / (1 + x(t - lag) / qm)")
print(" - Andy: x(t + dt) = r(t) + (x(t) - r) * exp(-dt * death), where r(t) = birth * x(t - lag) / death / (1 + x(t - lag) / qm)")
###  Make these floats, not integers!
death = 1.0 # 1/days
birth = 2.0 # 1/days
qm = 6.0 # population controller
dt = 1.0 # days
lag = 10.0 # days: best that this be a multiple of dt
x0 = 1 # population before t = 0
end_time = 10 * lag # days
print("Parameters produce steady-state x = " + str(qm * (birth / death - 1)) + "\n")

times = np.arange(-lag, end_time, dt)
x_euler = np.array([float(x0) for t in times])
x_andy = np.array([float(x0) for t in times])


def euler(x, t_plus_dt):
    t_plus_dt_index = np.argmin(np.abs(times - t_plus_dt))
    t_index = np.argmin(np.abs(times - (t_plus_dt - dt)))
    lagged_index = np.argmin(np.abs(times - (t_plus_dt - dt - lag)))
    return (t_plus_dt_index, x[t_index] + dt * (-death * x[t_index] + birth * x[lagged_index] / (1.0 + x[lagged_index] / qm)))

def andy(x, t_plus_dt):
    t_plus_dt_index = np.argmin(np.abs(times - t_plus_dt))
    t_index = np.argmin(np.abs(times - (t_plus_dt - dt)))
    lagged_index = np.argmin(np.abs(times - (t_plus_dt - dt - lag)))
    r = birth * x[lagged_index] / death / (1 + x[lagged_index] / qm)
    return (t_plus_dt_index, r + (x[t_index] - r) * np.exp(-dt * death))

for t in np.arange(dt, end_time, dt):
    ind, expl = euler(x_euler, t)
    x_euler[ind] = expl
    
for t in np.arange(dt, end_time, dt):
    ind, expl = andy(x_andy, t)
    x_andy[ind] = expl
    
plt.figure()
plt.plot(times, x_euler, label = "explicit Euler")
plt.plot(times, x_andy, label = "Andy")
plt.legend()
plt.grid()
plt.xlim([0, plt.xlim()[1]])
plt.xlabel("t")
plt.ylabel("x")
plt.title("Death=" + str(death) + " Birth=" + str(birth) + " qm=" + str(qm) + " lag=" + str(lag) + " dt=" + str(dt))
plt.show()
plt.close()

sys.exit(0)
