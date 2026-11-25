from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("Run without -O so that the rational assertions are active.")
sys.path.insert(0, str(Path(__file__).resolve().parent / "certificate"))
from verificar_p_general import verificar_perfil
from verificar_H5A_fusiones_sin_repetir import general

steps, sites = verificar_perfil()
print("Action and flipped-site bounds per data source:", steps, sites, flush=True)
general()
