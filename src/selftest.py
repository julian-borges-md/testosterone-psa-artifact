"""Proves verify.py can fail. Copies results to a temporary directory, perturbs
one committed estimate, and requires verify.py to exit non-zero."""
import shutil, subprocess, sys, tempfile, pandas as pd, os
tmp = tempfile.mkdtemp()
shutil.copytree('results', os.path.join(tmp, 'results'))
p = os.path.join(tmp, 'results/stage1/mr_results.tsv')
d = pd.read_csv(p, sep='\t'); d.loc[(d.analysis == 'biot->psa') & (d.method == 'IVW'), 'beta'] = 0.0
d.to_csv(p, sep='\t', index=False)
r = subprocess.run([sys.executable, 'src/verify.py', os.path.join(tmp, 'results/')], capture_output=True, text=True)
shutil.rmtree(tmp)
if r.returncode == 0:
    print('SELFTEST FAILED: verify.py did not detect a perturbed estimate', file=sys.stderr); sys.exit(1)
print('Selftest passed: verify.py detected the perturbation and exited', r.returncode)
