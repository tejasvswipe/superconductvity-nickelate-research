"""Run from repository root with: python finite_size_control.py (NumPy required)."""
import csv
from pathlib import Path
import numpy as np
from spatial_covariance_sensitivity import spectral_density, spans, wilson

ROOT = Path(__file__).resolve().parent
FAMILIES = ("gaussian_spectrum", "lorentzian_spectrum")
LAGS = (1, 2, 3, 4, 5, 6, 7, 8, 12, 16)


def kernel(n, family):
    m = 3*n
    lo, hi = 0., 8.
    for _ in range(48):
        ell = (lo+hi)/2
        power = spectral_density(family, ell, m)
        covariance = np.fft.irfft2(power[:, :m//2+1], s=(m,m))
        if covariance[0,1]/covariance[0,0] < .8:
            lo = ell
        else:
            hi = ell
    return np.sqrt(power[:, :m//2+1])


def run():
    out = ROOT/'finite_size_results'
    out.mkdir(exist_ok=True)
    with (out/'maps.csv').open('w',newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['size','family','batch','trial','active_sites','edge','edge_shuffle','diagonal','diagonal_shuffle']+[f'cov_lag_{d}' for d in LAGS])
        for n in (32,64,128):
            for fi, family in enumerate(FAMILIES):
                amplitude = kernel(n, family)
                for batch in (1,2):
                    field_rng, shuffle_rng = [np.random.default_rng(s) for s in np.random.SeedSequence([784203,n,fi,batch]).spawn(2)]
                    for trial in range(600):
                        noise = field_rng.standard_normal((3*n,3*n))
                        crop = np.fft.irfft2(np.fft.rfft2(noise)*amplitude,s=noise.shape)[n:2*n,n:2*n]
                        z = (crop-crop.mean())/crop.std()
                        active = abs(z) <= .08/.13
                        shuffled = shuffle_rng.permutation(active.ravel()).reshape(n,n)
                        assert active.sum()==shuffled.sum()
                        counts = [int(spans(mask,diagonal)) for diagonal in (False, True) for mask in (active,shuffled)]
                        assert counts[0]<=counts[2] and counts[1]<=counts[3]
                        cov = [(np.mean(z[:, :-d]*z[:, d:])+np.mean(z[:-d,:]*z[d:,:]))/2 for d in LAGS]
                        writer.writerow([n,family,batch,trial,active.sum(),*counts,*cov])
                    print(n,family,batch,flush=True)
    summarize(out)


def summarize(out):
    with (out/'maps.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    with (out/'summary.csv').open('w',newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['size','family','batch','maps','edge','edge_shuffle','diagonal','diagonal_shuffle','edge_ci_low','edge_ci_high','edge_only_original','edge_only_shuffle']+[f'mean_cov_lag_{d}' for d in LAGS])
        for n in (32,64,128):
            for family in FAMILIES:
                for batch in ('1','2','pooled'):
                    group=[r for r in rows if int(r['size'])==n and r['family']==family and (batch=='pooled' or batch==r['batch'])]
                    counts=[sum(int(r[k]) for r in group) for k in ('edge','edge_shuffle','diagonal','diagonal_shuffle')]
                    pairs=[sum(r['edge']==a and r['edge_shuffle']==b for r in group) for a,b in [('1','0'),('0','1')]]
                    curve=[np.mean([float(r[f'cov_lag_{d}']) for r in group]) for d in LAGS]
                    writer.writerow([n,family,batch,len(group),*counts,*wilson(counts[0],len(group)),*pairs,*curve])


if __name__=='__main__':
    run()
