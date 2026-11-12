"""Diseño y evaluación reproducible de seis filtros causales FIR e IIR.

Uso: python filtros_digitales.py --salida resultados
Dependencias: NumPy, SciPy y Matplotlib. No requiere archivos externos.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FS = 2000
DURATION = 4.0
SEED = 2026
SIGMA = 0.5
TONES = (50, 200, 400)
AMPLITUDES = (1.0, 0.8, 0.6)
CASES = {
    'pasa_bajos': {'title': 'Pasa bajos', 'cutoff': 100, 'btype': 'lowpass', 'pass_zero': True, 'target': 0},
    'pasa_altos': {'title': 'Pasa altos', 'cutoff': 300, 'btype': 'highpass', 'pass_zero': False, 'target': 2},
    'pasa_bandas': {'title': 'Pasa bandas', 'cutoff': [150, 250], 'btype': 'bandpass', 'pass_zero': False, 'target': 1},
}
COLORS = {'FIR': '#137c75', 'IIR': '#d46628'}


def design(kind, family, size):
    """FIR: size coeficientes Hamming. IIR: size orden de prototipo Butterworth."""
    c = CASES[kind]
    if family == 'FIR':
        return signal.firwin(size, c['cutoff'], fs=FS, window='hamming', pass_zero=c['pass_zero'])
    return signal.butter(size, c['cutoff'], btype=c['btype'], fs=FS, output='sos')


def apply(coeff, family, x):
    """Una pasada causal, condiciones iniciales nulas; no se usa filtfilt."""
    return signal.lfilter(coeff, [1.0], x) if family == 'FIR' else signal.sosfilt(coeff, x)


def response(coeff, family, freq):
    if family == 'FIR':
        return signal.freqz(coeff, worN=freq, fs=FS)[1]
    return signal.sosfreqz(coeff, worN=freq, fs=FS)[1]


def db(x):
    """Magnitud a dB; límite solo para evitar log(0) en gráficas."""
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12))


def power(x):
    return float(np.mean(np.asarray(x)**2))


def snr(desired, unwanted):
    return float(10 * np.log10(power(desired)/power(unwanted)))


def masks(kind, f):
    # Bandas de evaluación independientes de las frecuencias de corte.
    if kind == 'pasa_bajos':
        return f <= 70, f >= 150
    if kind == 'pasa_altos':
        return f >= 400, f <= 200
    return (f >= 185) & (f <= 215), (f <= 100) | (f >= 300)


def measure(kind, family, size, coeff, t, components, noise):
    target = CASES[kind]['target']
    desired = components[target]
    unwanted = sum(x for i, x in enumerate(components) if i != target) + noise
    yd, yu = apply(coeff, family, desired), apply(coeff, family, unwanted)
    y = apply(coeff, family, desired + unwanted)
    # Descartar 0.5 s de arranque. Se mantiene el retardo físico en las señales.
    valid = t >= 0.5
    f = np.linspace(0, FS/2, 20001)
    h = response(coeff, family, f)
    pb, sb = masks(kind, f)
    gain = db(h)
    tone_db = db(response(coeff, family, np.array(TONES)))
    # Retardo de grupo por derivada de fase respecto a frecuencia angular.
    gd = -np.gradient(np.unwrap(np.angle(h)), 2*np.pi*f/FS)
    idx = np.argmin(abs(f-TONES[target]))
    poles = np.array([]) if family == 'FIR' else signal.sos2zpk(coeff)[1]
    row = dict(tipo=kind, familia=family, parametro=size,
               orden_real=size-1 if family=='FIR' else size*(2 if kind=='pasa_bandas' else 1),
               snr_entrada_db=snr(desired[valid], unwanted[valid]),
               snr_salida_db=snr(yd[valid], yu[valid]),
               perdida_paso_max_db=float(-min(gain[pb])),
               rechazo_min_db=float(-max(gain[sb])),
               ganancia_50_db=float(tone_db[0]), ganancia_200_db=float(tone_db[1]),
               ganancia_400_db=float(tone_db[2]), retardo_objetivo_ms=float(gd[idx]/FS*1000),
               radio_polo_max=float(max(abs(poles))) if len(poles) else 0.0)
    row['mejora_snr_db'] = row['snr_salida_db']-row['snr_entrada_db']
    row['cumple_1db_20db'] = row['perdida_paso_max_db'] <= 1 and row['rechazo_min_db'] >= 20
    assert np.all(np.isfinite(y))
    assert np.max(abs(y-yd-yu)) < 1e-10, 'Fallo de superposición'
    assert row['radio_polo_max'] < 1, 'Filtro IIR inestable'
    return row, y, yd


def savefig(fig, path):
    fig.savefig(path, dpi=180, bbox_inches='tight'); plt.close(fig)


def run(out: Path):
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    t = np.arange(int(FS*DURATION))/FS
    components = [a*np.sin(2*np.pi*f*t) for f, a in zip(TONES, AMPLITUDES)]
    noise = np.random.default_rng(SEED).normal(0, SIGMA, t.size)
    clean = sum(components); x = clean+noise
    rows, baseline, coeffs = [], {}, {}
    for kind in CASES:
        baseline[kind] = {}
        for family, sizes in [('FIR', [51, 101, 201]), ('IIR', [2, 4, 6])]:
            for size in sizes:
                coeff = design(kind, family, size)
                row, y, yd = measure(kind, family, size, coeff, t, components, noise)
                rows.append(row)
                if size == (201 if family=='FIR' else 4):
                    baseline[kind][family] = {'metrics': row, 'y': y, 'desired': yd}
                    coeffs[kind+'_'+family] = coeff.tolist()
    with (out/'metricas_todas_configuraciones.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    base_rows=[baseline[k][fam]['metrics'] for k in CASES for fam in ['FIR','IIR']]
    with (out/'metricas_principales.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(base_rows)
    (out/'coeficientes.json').write_text(json.dumps(coeffs,indent=2),encoding='utf8')
    versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)
    (out/'entorno.json').write_text(json.dumps(versions,indent=2),encoding='utf8')
    valid = t >= .5
    f, psdx = signal.welch(x[valid],fs=FS,window='hann',nperseg=2000,noverlap=1000,scaling='density')
    fig, ax = plt.subplots(2,1,figsize=(10,6),layout='constrained')
    sel=(t>=.5)&(t<=.6)
    ax[0].plot(t[sel],x[sel],color='#8a8a8a',label='Entrada con ruido',lw=1)
    ax[0].plot(t[sel],clean[sel],color='#264e70',label='Suma de tonos sin ruido',lw=1.3)
    ax[0].set(xlabel='Tiempo (s)',ylabel='Amplitud (u.a.)',title='Señal de entrada'); ax[0].legend(loc='upper right',fontsize=9)
    ax[1].plot(f,10*np.log10(np.maximum(psdx,1e-15)),color='#264e70')
    ax[1].set(xlabel='Frecuencia (Hz)',ylabel='PSD (dB u.a.²/Hz)',xlim=(0,1000),ylim=(-80,10),title='Densidad espectral de potencia por Welch')
    for fr in TONES: ax[1].axvline(fr,color='#aaaaaa',ls=':',lw=.8)
    savefig(fig,out/'01_entrada.png')
    grid=np.linspace(0,1000,5001)
    fig, axs=plt.subplots(3,1,figsize=(10,8),layout='constrained')
    for ax,kind in zip(axs,CASES):
        for fam,size in [('FIR',201),('IIR',4)]: ax.plot(grid,db(response(design(kind,fam,size),fam,grid)),label=fam,color=COLORS[fam])
        ax.set(xlim=(0,1000),ylim=(-100,5),ylabel='Ganancia (dB)',title=CASES[kind]['title']); ax.grid(alpha=.2); ax.legend(loc='lower right')
    axs[-1].set_xlabel('Frecuencia (Hz)'); savefig(fig,out/'02_respuestas.png')
    for ix,kind in enumerate(CASES,3):
        fig,ax=plt.subplots(2,1,figsize=(10,6.3),layout='constrained')
        target=components[CASES[kind]['target']]
        ax[0].plot(t[sel],x[sel],color='#bbbbbb',lw=.7,label='Entrada')
        for fam in ['FIR','IIR']: ax[0].plot(t[sel],baseline[kind][fam]['y'][sel],color=COLORS[fam],label=f'{fam} causal',lw=1.2)
        # La referencia sin filtrar hace visible el desplazamiento de fase.
        ax[0].plot(t[sel],target[sel],color='black',ls=':',lw=1,label='Tono objetivo sin filtrar')
        ax[0].set(xlabel='Tiempo (s)',ylabel='Amplitud (u.a.)',title=CASES[kind]['title']+' en el tiempo'); ax[0].legend(loc='upper right',fontsize=8,ncols=2)
        ax[1].plot(f,10*np.log10(np.maximum(psdx,1e-15)),color='#aaaaaa',label='Entrada')
        for fam in ['FIR','IIR']:
            fy,py=signal.welch(baseline[kind][fam]['y'][valid],fs=FS,window='hann',nperseg=2000,noverlap=1000,scaling='density')
            ax[1].plot(fy,10*np.log10(np.maximum(py,1e-15)),color=COLORS[fam],label=fam)
        ax[1].set(xlabel='Frecuencia (Hz)',ylabel='PSD (dB u.a.²/Hz)',xlim=(0,600),ylim=(-100,10),title='Entrada y salidas en frecuencia'); ax[1].legend()
        savefig(fig,out/f'{ix:02d}_{kind}.png')
    fig,ax=plt.subplots(figsize=(10,4),layout='constrained'); pos=np.arange(3)
    for j,fam in enumerate(['FIR','IIR']): ax.bar(pos+(j-.5)*.32,[baseline[k][fam]['metrics']['snr_salida_db'] for k in CASES],width=.32,label=fam,color=COLORS[fam])
    ax.set_xticks(pos,[CASES[k]['title'] for k in CASES]); ax.set_ylabel('SNR de salida (dB)'); ax.legend(); ax.grid(axis='y',alpha=.2); savefig(fig,out/'06_snr.png')
    fig,axs=plt.subplots(3,1,figsize=(10,7),layout='constrained')
    for ax,kind in zip(axs,CASES):
        pb,_=masks(kind,grid)
        for fam,size in [('FIR',201),('IIR',4)]:
            hh=response(design(kind,fam,size),fam,grid); gd=-np.gradient(np.unwrap(np.angle(hh)),2*np.pi*grid/FS)/FS*1000
            ax.plot(grid[pb],gd[pb],label=fam,color=COLORS[fam])
        ax.set(title=CASES[kind]['title'],ylabel='Retardo (ms)'); ax.legend(); ax.grid(alpha=.2)
    axs[-1].set_xlabel('Frecuencia en la banda de paso (Hz)'); savefig(fig,out/'07_retardo.png')
    # Datos ligeros para gráficos editables en la presentación.
    fg=np.linspace(0,600,301)
    deck={'fs':FS,'duration':DURATION,'seed':SEED,'versions':versions,'rows':base_rows,'all_rows':rows,'freq':fg.tolist(),
          'responses':{k:{fam:np.maximum(db(response(design(k,fam,size),fam,fg)),-100).tolist() for fam,size in [('FIR',201),('IIR',4)]} for k in CASES},
          'time':t[sel].tolist(),'input':x[sel].tolist(),'clean':clean[sel].tolist(),
          'time_bp':{fam:baseline['pasa_bandas'][fam]['y'][sel].tolist() for fam in ['FIR','IIR']}}
    (out/'resumen.json').write_text(json.dumps(deck,ensure_ascii=False,indent=2),encoding='utf8')
    print('Simulación completa: 18 configuraciones, 6 filtros principales y 7 gráficas.')
    for row in base_rows: print(row['tipo'],row['familia'],'SNR:',round(row['snr_salida_db'],2),'rechazo:',round(row['rechazo_min_db'],2),'cumple:',row['cumple_1db_20db'])
    return deck

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--salida',type=Path,default=Path('resultados'))
    run(parser.parse_args().salida)
