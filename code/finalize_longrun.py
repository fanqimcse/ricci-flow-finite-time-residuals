"""Wait for every ricci_flow_longrun worker to finish, then refresh the
analysis tables, figures and the HTML report.

Run it once and walk away:

    python finalize_longrun.py

It polls the process table every 30 seconds, and when no worker is left it
runs ``analyze_longrun.py`` followed by ``make_longrun_report.py``.
"""

from __future__ import print_function

import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
POLL_SECONDS = 30
MAX_WAIT_SECONDS = 12 * 3600


PS = ("Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | "
      "Where-Object { $_.CommandLine -like '*ricci_flow_longrun*' } | "
      "Measure-Object | Select-Object -ExpandProperty Count")


def count_workers():
    """Number of longrun workers, matched by command line (not just by name).

    Counting by image name is wrong: the shell that launches a worker and the
    watcher itself are also python.exe / bash.exe processes.
    """
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-Command", PS],
            capture_output=True, text=True, timeout=120
        ).stdout.strip()
    except Exception:
        return -1
    try:
        return int(out.splitlines()[-1])
    except (ValueError, IndexError):
        return -1


def main():
    started = time.time()
    log = open(os.path.join(BASE, "_finalize_log.txt"), "w", encoding="utf-8")

    def say(message):
        line = "[{}] {}".format(time.strftime("%H:%M:%S"), message)
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    say("watcher started, waiting for longrun workers to exit")
    while time.time() - started < MAX_WAIT_SECONDS:
        alive = count_workers()
        if alive == 0:
            say("all workers finished")
            break
        if alive < 0:
            say("could not read the process table; assuming workers are done")
            break
        if (time.time() - started) % 600 < POLL_SECONDS:
            say("still waiting (python.exe processes seen: {})".format(alive))
        time.sleep(POLL_SECONDS)
    else:
        say("watch timeout reached; finalising with whatever is on disk")

    for script in ("analyze_longrun.py", "make_longrun_report.py"):
        path = os.path.join(BASE, script)
        say("running {}".format(script))
        result = subprocess.run([PY, "-u", path], cwd=BASE,
                                capture_output=True, text=True)
        say("{} exit={}".format(script, result.returncode))
        if result.returncode != 0:
            say(result.stdout[-2000:])
            say(result.stderr[-2000:])

    say("DONE")


if __name__ == "__main__":
    main()
