import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_left
import time
import pdb

__all__ =  ['lookup',
            'triangular_hyetograph',
            'rectangular_hyetograph',
            'time_area_hydrograph',
            'calculate_storage_volumes',
            'balance',
            'compute_calcs'
            ]

def lookup(x, xs, ys, interpolate=True):
    """Interpolate"""

    if x <= xs[0]:  return ys[0]
    if x >= xs[-1]: return ys[-1]

    i = bisect_left(xs, x)
    if interpolate:
        k = (x - xs[i-1])/(xs[i] - xs[i-1])
        y = k*(ys[i]-ys[i-1]) + ys[i-1]
    else:
        y = ys[i-1]

    return y

def triangular_hyetograph(t, R, D):
    """Triangular-shaped hyetograph"""

    d_peak = D/2
    i_peak = R/d_peak
    x = np.zeros(len(t))
    left_mask = np.where(t <= D/2)
    right_mask = np.where(t > D/2)
    x[left_mask] = i_peak * (t[left_mask]/d_peak)
    x[right_mask] = np.where((1 - (t[right_mask]-d_peak)/d_peak) > 0, i_peak * (1 - (t[right_mask]-d_peak)/d_peak), 0)
    return x

def rectangular_hyetograph(t, R, D):
    """Rectangular-shaped hyetograph"""

    x = np.zeros(len(t))
    mask = np.where(t <= D)
    x[mask] = R/D
    return x

def time_area_hydrograph(x, tc, w=[0.25, 0.25, 0.25, 0.25], tstep=1):
    """Time area hydrograph"""
    tc = tc / tstep
    z = np.zeros(len(x) + int(max(tc)))
    for i in range(len(tc)):
        z[int(tc[i]):int(tc[i])+len(x)] += x*w[i]
    return z

def calculate_storage_volumes(l_pond, a_pond):
    """Calculates storage volumes"""
    v_pond = []
    for i in range(len(a_pond)):
        v = np.trapz(a_pond[:i+1], x=l_pond[:i+1])
        v_pond.append(v)
    return v_pond

def balance(inputs, init=0, storage_area=2000, area=None, rc=1, 
            constant_outflow=False, outflow_rate=10, outflow_table=None, 
            interpolate=True, tstep = 1, l_pond=None, v_pond=None):

    """Balance"""
    
    s = np.zeros(len(inputs))
    l = np.zeros(len(inputs))
    o = np.zeros(len(inputs))
    i = np.zeros(len(inputs))
    
    for j in range(len(s)):
        
        i[j] = inputs[j] * area * rc
        
        if j>0:
            if constant_outflow:
                o[j] = outflow_rate
            else:
                o[j] = lookup(l[j-1], outflow_table['x'], outflow_table['y'], interpolate=interpolate)
            
            s[j] = max(0, s[j-1] + (i[j] - o[j]) * tstep / 1000)
        
        else:
            s[j] = max(0, init + (i[j] - o[j]) * tstep / 1000)
        
        l[j] = lookup(x=s[j], xs=v_pond, ys=l_pond)
            
    return s, l, i, o

def compute_calcs(df, RPs, area, rc, pump_rate, constant_outflow=False, 
                  outflow_table=None, interpolate=True, storage_area=2000, 
                  tstep=1, tc=None, w=None, cv_winter=0.84, l_pond=None, 
                  a_pond=None, v_pond=None, calculate_volume=True):

    """Compute calculations"""
    
    durations = df.index.values
    results = {}
    design_volumes, design_times, design_levels, design_outflows, design_durations = [], [], [], [], []
    
    if calculate_volume:
        v_pond = calculate_storage_volumes(l_pond, a_pond)
    else:
        pass
    
    for RP in RPs[:]:
        results[f'{RP} year'] = {}
        peak_storage_volumes, peak_storage_times, peak_storage_levels, peak_storage_outflows, peak_storage_durations = [], [], [], [], []
        for i, duration in enumerate(durations):
            
            results[f'{RP} year'][f'{duration} min'] = {}
            d = results[f'{RP} year'][f'{duration} min']
            d['R'] = df[f'{RP} year rainfall (mm)'].loc[duration] * (1 + cc_factor) * cv_winter
            d['D'] = duration*60
                        
            d['time'] = np.arange(0, 2*d['D']+max(tc), tstep)
            
            d['intensity'] = triangular_hyetograph(t=d['time'], R=d['R'], D=d['D'])
            #d['intensity'] = rectangular_hyetograph(t=d['time'], R=d['R'], D=d['D'])
            
            d['hydrograph'] = time_area_hydrograph(d['intensity'], tc=tc, w=w, tstep=tstep)
            
            d['storage'], d['level'], d['inflow'], d['outflow'] = balance(inputs=d['hydrograph'], 
                                                                          outflow_rate=pump_rate, 
                                                                          area=area, 
                                                                          rc=rc, 
                                                                          constant_outflow=constant_outflow, 
                                                                          outflow_table=outflow_table,
                                                                          interpolate=interpolate,
                                                                          storage_area=storage_area,
                                                                          tstep=tstep,
                                                                          l_pond=l_pond,
                                                                          v_pond=v_pond)
            

            d['peak_storage_volume'] = np.max(d['storage'])
            peak_storage_volumes.append(d['peak_storage_volume'])
            
            peak_storage_index = np.where(d['storage']==d['peak_storage_volume'])[0]
            
            d['peak_storage_time'] = int(d['time'][peak_storage_index])
            peak_storage_times.append(d['peak_storage_time'])
            
            d['peak_storage_level'] = float(d['level'][peak_storage_index])
            peak_storage_levels.append(d['peak_storage_level'])
            
            d['peak_storage_outflow'] = float(d['outflow'][peak_storage_index])
            peak_storage_outflows.append(d['peak_storage_outflow'])
            
            peak_storage_durations.append(duration)

        #results[f'{RP} year'][f'{duration} min'] = d
        
        results[f'{RP} year']['peak_storage_volumes'] = peak_storage_volumes
        results[f'{RP} year']['peak_storage_times'] = peak_storage_times
        results[f'{RP} year']['peak_storage_levels'] = peak_storage_levels
        results[f'{RP} year']['peak_storage_outflows'] = peak_storage_outflows

        design_index = int(np.where(peak_storage_volumes==max(peak_storage_volumes))[0])
        
        design_volumes.append(peak_storage_volumes[design_index])
        design_times.append(peak_storage_times[design_index])
        design_levels.append(peak_storage_levels[design_index])
        design_outflows.append(peak_storage_outflows[design_index])
        design_durations.append(peak_storage_durations[design_index])
        
        #design_level = int(np.where(peak_storage_volumes==max(peak_storage_volumes))[0])


    results['design_volumes'] = design_volumes
    results['design_times'] = design_times
    results['design_levels'] = design_levels
    results['design_outflows'] = design_outflows
    results['design_durations'] = design_durations
    
    return results