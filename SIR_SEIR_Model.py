# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 18:40:51 2020

@author: Lissette 
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Población total, N
N =   17e6
# Valores iniciales de personas Infectadas (I0) y personas que se recuperaron (R0)
I0, R0 = 1, 0
# El resto de la población, S0 son las personas que están sujetos a infección inicialmente.
S0 = N - I0 - R0
# Tasa de transmisión, tasa de recuperación en 1/días
#beta= R0 * gamma = 1.46 * 1/14 = 0.1042
#gamma= 1/14
#R0= 1.46 (tranmition in Ecuador)
beta, gamma = 0.1042, 1/14

# Puntos en la gráfica (En días)
t = np.linspace(0, 1000, 500)

# Las ecuaciones diferenciales del modelo SIR
def deriv(y, t, N, beta, gamma):
   S, I, R = y
   dSdt = -beta * S * I / N
   dIdt = beta * S * I / N - gamma * I
   dRdt = gamma * I
   return dSdt, dIdt, dRdt

# Vector de las condiciones iniciales
y0 = S0, I0, R0
# Resolver el sistema de ecuaciones diferenciales, en la secuencia de días que ya definimos
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/17e6, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I/17e6, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/17e6, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number of people (million)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()


#-----------------------ON GITHUB DINAMYC GRAPH----------------------------
import plotly.graph_objects as go

fig = go.Figure()
fig.add_trace(go.Scatter(x=t, y=I,mode='lines',name='Infected'))
fig.show()

#-----------------------ON GITHUB DINAMYC GRAPH----------------------------

#parametros en varias graph
import numpy as np
from scipy.integrate import odeint

valores = [(2.2,1/2.3), (5.1,1/7), (1,1/3)]

def modelo(beta, gamma):
 # Población total, N
      N = 5.2e6
 # Valores iniciales de personas Infectadas (I0) y personas que se recuperaron (R0)
      I0, R0 = 1, 0
 # El resto de la población, S0 son las personas que están sujetos a infección inicialmente.
      S0 = N - I0 - R0
 # Puntos en la gráfica (En días)
      t = np.linspace(0, 40, 40)

# Las ecuaciones diferenciales del modelo SIR
def deriv(y, t, N, beta, gamma):
      S, I, R = y
      dSdt = -beta * S * I / N
      dIdt = beta * S * I / N - gamma * I
      dRdt = gamma * I
      return dSdt, dIdt, dRdt
# Vector de las condiciones iniciales
      y0 = S0, I0, R0
  # Resolver el sistema de ecuaciones diferenciales, en la secuencia de días que ya definimos
      ret = odeint(deriv, y0, t, args=(N, beta, gamma))
      S, I, R = ret.T
      return (I)

import cufflinks as cf
import plotly.offline as py

#ciclo for para cada valor en valores
py.iplot([{
   'x': t,
   'y': modelo(*valor),
   'name': str(valor),
}  for valor in valores], filename='cufflinks/multiple-lines-on-same-chart')

    
#Reference
#https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/