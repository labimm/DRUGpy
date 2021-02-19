import subprocess
import sys


def pip(args):
    process = subprocess.Popen(
        [sys.executable, "-m", "pip", "--disable-pip-version-check"] + args,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = process.communicate()
    return process, out, err


try:
    from pymol_labimm import init_plugin
except:
    proc, out, err = pip(["install", "--upgrade", "pymol-labimm"])
    if out:
        print(out)
    if err:
        print(err)
    from pymol_labimm import init_plugin


init_plugin()
