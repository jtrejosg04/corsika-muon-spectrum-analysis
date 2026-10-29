#!/usr/bin/env python3
"""Compute the atmospheric muon spectrum using Wiebel-Sooth.

Run with:
    python3 compute_wiebel_sooth_muon_spectrum.py

Input files must follow the directory structure
documented in README.md.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
RUN = BASE.parent
KEYS = ['run_id', 'shower_id']
GAMMA_GENERACION = 2.7
OMEGA_PRIM = np.pi * np.sin(np.deg2rad(30.0))**2
OMEGA_MU = np.pi * (1.0 - 0.975**2)
# (símbolo, Z, A, phi0, gamma), de prueba_B6_1_nucleones.py
NUCLEOS = [
    ('H',1,1,11.51e-2,2.77), ('He',2,4,7.19e-2,2.64),
    ('Li',3,7,2.08e-3,2.54), ('Be',4,9,4.74e-4,2.75),
    ('B',5,11,8.95e-4,2.95), ('C',6,12,1.06e-2,2.66),
    ('N',7,14,2.35e-3,2.72), ('O',8,16,1.57e-2,2.68),
    ('F',9,19,3.28e-4,2.69), ('Ne',10,20,4.60e-3,2.64),
    ('Na',11,23,7.54e-4,2.66), ('Mg',12,24,8.01e-3,2.64),
    ('Al',13,27,1.15e-3,2.66), ('Si',14,28,7.96e-3,2.75),
    ('P',15,31,2.70e-4,2.69), ('S',16,32,2.29e-3,2.55),
    ('Cl',17,35,2.94e-4,2.68), ('Ar',18,40,8.36e-4,2.64),
    ('K',19,39,5.36e-4,2.65), ('Ca',20,40,1.47e-3,2.70),
    ('Sc',21,45,3.04e-4,2.64), ('Ti',22,48,1.13e-3,2.61),
    ('V',23,51,6.31e-4,2.63), ('Cr',24,52,1.36e-3,2.67),
    ('Mn',25,55,1.35e-4,2.46), ('Fe',26,56,1.78e-2,2.60),
    ('Co',27,59,7.51e-5,2.72), ('Ni',28,58,9.96e-4,2.51),
]
WS = {s: (z, a, phi, gam) for s, z, a, phi, gam in NUCLEOS}
# (runs, N_total_banda, E_min, E_max, A_simulada)
PROTONES = [
    ([201,202,203],50500,10,100,1),
    ([204,205],10500,100,1000,1),
    ([206,207],5200,1000,10000,1),
    ([208,209],550,10000,100000,1),
]
HELIO = [
    ([210],10000,10,100,4),
    ([211],10000,100,1000,4),
    ([212],5000,1000,10000,4),
    ([213],500,10000,60000,4),
]

def j_ws(E_nucleo, simbolo):
    """Intensidad por GeV de energía total nuclear; cero fuera de E>=10Z."""
    z, a, phi, gamma = WS[simbolo]
    E = np.asarray(E_nucleo, dtype=float)
    J = np.zeros_like(E)
    valido = E >= 10.0*z
    J[valido] = phi/1000.0 * (E[valido]/1000.0)**(-gamma)
    return J

def densidad_generacion(E, emin, emax):
    g = GAMMA_GENERACION
    C = (1.0-g)/(emax**(1.0-g)-emin**(1.0-g))
    return C*np.asarray(E)**(-g)

def cargar_grupo(runs, n_esperado, emin, emax, especie):
    prims, mus = [], []
    for run in runs:
        p = pd.read_csv(BASE / 'catalogos' / f'primarios_{run}.csv')
        m = pd.read_csv(BASE / 'catalogos' / f'muones_{run}.csv')
        assert not p.duplicated(['shower_id']).any(), f'Primarios duplicados run {run}'
        p['run_id'] = run
        m['run_id'] = run
        prims.append(p)
        mus.append(m)
    p = pd.concat(prims, ignore_index=True)
    m = pd.concat(mus, ignore_index=True)
    if len(p) != n_esperado:
        raise ValueError(f'{especie} runs {runs}: N={len(p)} vs esperado {n_esperado}')
    E = p['energy_GeV'].to_numpy(dtype=float)
    if not np.all(np.isfinite(E) & (E>0)):
        raise ValueError('Energías primarias inválidas')
    if E.min() < emin*(1-1e-8) or E.max() > emax*(1+1e-8):
        raise ValueError(f'Energías fuera de banda en runs {runs}')
    g = densidad_generacion(E, emin, emax)
    if especie == 'protones':
        Jp = j_ws(E, 'H')
        Jheavy = np.zeros_like(E)
        for simbolo, z, a, phi, gamma in NUCLEOS:
            if z >= 3:
                Jheavy += a*a*j_ws(a*E, simbolo)
        p['w_p'] = OMEGA_PRIM*Jp/(n_esperado*g)
        p['w_he'] = 0.0
        p['w_pesados'] = OMEGA_PRIM*Jheavy/(n_esperado*g)
    else:
        p['w_p'] = 0.0
        p['w_he'] = OMEGA_PRIM*j_ws(E, 'He')/(n_esperado*g)
        p['w_pesados'] = 0.0
    for col in ['w_p','w_he','w_pesados']:
        if not np.all(np.isfinite(p[col]) & (p[col]>=0)):
            raise ValueError(f'Pesos inválidos: {col}')
    # Solo datos de muones y multiplicidad, sin usar pesos previos
    needed = ['run_id','shower_id','nivel','p','cos_theta']
    m = m[needed].copy()
    if not pd.MultiIndex.from_frame(m[KEYS]).isin(pd.MultiIndex.from_frame(p[KEYS])).all():
        raise ValueError(f'Muones sin primario en {runs}')
    return p, m

def estadistica_grupo(p, m, lo, hi):
    aceptados = m[(m['nivel']==1) & (m['p']>=lo) & (m['p']<hi) &
                  (m['cos_theta']>0.975) & (m['cos_theta']<1.0)]
    counts = aceptados.groupby(KEYS).size()
    indice = pd.MultiIndex.from_frame(p[KEYS])
    n = counts.reindex(indice, fill_value=0).to_numpy(dtype=float)
    yp = p['w_p'].to_numpy()*n
    yh = p['w_he'].to_numpy()*n
    ypes = p['w_pesados'].to_numpy()*n
    ytotal = yp + yh + ypes  # conserva covarianza protones-pesados
    def sv(y):
        N = len(y)
        if N <= 1:
            raise ValueError('Grupo con menos de dos lluvias')
        return float(y.sum()), float(N/(N-1)*np.sum((y-y.mean())**2))
    return {k:sv(y) for k,y in [('p',yp),('he',yh),('pesados',ypes),('total',ytotal)]}, int(n.sum())

def main():
    exp = pd.read_csv(RUN / 'L3C_datos_experimentales.csv')
    required = ['p_min_GeV','p_max_GeV','p_representativo_GeV']
    if not set(required).issubset(exp.columns):
        raise ValueError(f'Faltan columnas en L3C: {set(required)-set(exp.columns)}')
    grupos = ([cargar_grupo(args[0], args[1], args[2], args[3], 'protones') for args in PROTONES] +
              [cargar_grupo(args[0], args[1], args[2], args[3], 'helio') for args in HELIO])
    resultados = []
    for r in exp.itertuples(index=False):
        lo, hi = float(r.p_min_GeV), float(r.p_max_GeV)
        if not (0 < lo < hi):
            raise ValueError(f'Intervalo inválido {lo}, {hi}')
        S = {k:0.0 for k in ['p','he','pesados','total']}
        V = {k:0.0 for k in S}
        nmu = 0
        for p,m in grupos:
            sv, n = estadistica_grupo(p,m,lo,hi)
            nmu += n
            for k,(s,v) in sv.items():
                S[k] += s
                V[k] += v
        norma = (hi-lo)*OMEGA_MU
        fila = {'p_min_GeV':lo,'p_max_GeV':hi,
                'p_GeV':float(r.p_representativo_GeV),'n_muones':nmu}
        for k in S:
            fila[f'Phi_{k}'] = S[k]/norma
            fila[f'sigma_{k}'] = np.sqrt(V[k])/norma
        resultados.append(fila)
    out = pd.DataFrame(resultados)
    if not np.allclose(out['Phi_total'],out[['Phi_p','Phi_he','Phi_pesados']].sum(axis=1)):
        raise AssertionError('La suma de componentes no coincide')
    csv = BASE / 'L3C_ATMOD7_WS_puro.csv'
    png = BASE / 'L3C_ATMOD7_WS_puro.png'
    for path in (csv,png):
        if path.exists():
            raise FileExistsError(f'No se sobrescribe {path}')
    # Intensidad experimental: solo si hay columna identificable; sin inventar unidades
    candidatos = ['Phi_experimental','phi_experimental','Phi','phi','intensidad','intensidad_experimental']
    col_exp = next((c for c in candidatos if c in exp.columns),None)
    # La comparación histórica guarda la intensidad experimental en este CSV.
    y_exp_fallback = None
    if col_exp is None:
        archivo_comparacion = RUN / 'comparacion_p_He_L3C.csv'
        if archivo_comparacion.exists():
            comp = pd.read_csv(archivo_comparacion)
            if 'Phi_experimental' in comp.columns and len(comp)==len(out):
                if 'p_GeV' not in comp.columns or np.allclose(comp['p_GeV'],out['p_GeV']):
                    y_exp_fallback = pd.to_numeric(comp['Phi_experimental'],errors='coerce').to_numpy()
    fig, ax = plt.subplots(figsize=(9,6))
    x = out['p_GeV'].to_numpy()
    if col_exp:
        y = pd.to_numeric(exp[col_exp],errors='coerce').to_numpy()
        mask = np.isfinite(y)&(y>0)
        ax.plot(x[mask],y[mask],'ko',label=f'L3+C experimental ({col_exp})',ms=4)
    elif y_exp_fallback is not None:
        mask = np.isfinite(y_exp_fallback)&(y_exp_fallback>0)
        ax.plot(x[mask],y_exp_fallback[mask],'ko',label='L3+C experimental',ms=4)
    else:
        print('AVISO: no se detectó intensidad experimental. Revisar los CSV de L3+C.')
        print('Columnas disponibles:',list(exp.columns))
    for key,label,fmt in [('p','Protones WS','-'),('he','Helio WS','--'),
                          ('pesados','Pesados WS (superposición)',':')]:
        y = out[f'Phi_{key}'].to_numpy()
        mask = np.isfinite(y)&(y>0)
        ax.plot(x[mask],y[mask],fmt,label=label)
    yt = out['Phi_total'].to_numpy()
    st = out['sigma_total'].to_numpy()
    mask = np.isfinite(yt)&(yt>0)&np.isfinite(st)
    ax.errorbar(x[mask],yt[mask],yerr=st[mask],fmt='s-',capsize=2,ms=4,
                label='Total WS (incertidumbre MC)')
    ax.set(xscale='log',yscale='log',xlabel='Momento del muón (GeV/c)',
           ylabel=r'Intensidad diferencial [$m^{-2}\,s^{-1}\,sr^{-1}\,(GeV/c)^{-1}$]',
           title='L3+C vs CORSIKA–UrQMD · ATMOD7 · Wiebel-Sooth')
    ax.grid(True,which='both',alpha=.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out.to_csv(csv,index=False)
    fig.savefig(png,dpi=220)
    plt.close(fig)
    print('Resultados:',csv,'\nGráfica:',png)
    print('Nota: pesados = superposición; corte WS E_A >= 10 Z GeV.')
    print('Nota: confirmar unidades del CSV experimental y convención de energy_GeV.')

if __name__ == '__main__':
    main()
