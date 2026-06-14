import os
import sys
from pathlib import Path

from streamlit.web import cli


def _bundle_root():
    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))


def _configure_gemelli_worker():
    if os.environ.get("GEMELLI_WORKER_PYTHON"):
        return

    exe_dir = Path(sys.executable).resolve().parent
    candidates = [
        exe_dir / "gemelli-standalone" / "python.exe",
        exe_dir / "gemelli-standalone" / "Scripts" / "python.exe",
        exe_dir / "gemelli-standalone" / "bin" / "python",
        _bundle_root() / "gemelli-standalone" / "python.exe",
        _bundle_root() / "gemelli-standalone" / "Scripts" / "python.exe",
        _bundle_root() / "gemelli-standalone" / "bin" / "python",
    ]
    for worker_python in candidates:
        if worker_python.exists():
            os.environ["GEMELLI_WORKER_PYTHON"] = str(worker_python)
            return


if __name__ == "__main__":
    app_root = _bundle_root()
    os.chdir(app_root)
    _configure_gemelli_worker()
    cli._main_run_clExplicit(
        file=str(app_root / "Home.py"),
        command_line="streamlit run",
        args=[],
    )
