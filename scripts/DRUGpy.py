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
    from drugpy import init_plugin
except:
    proc, out, err = pip(["install", "--upgrade", "DRUGpy"])
    if out:
        print(out)
    if err:
        print(err)
    from drugpy import init_plugin


init_plugin()
