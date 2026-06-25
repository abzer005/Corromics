import os
import socket
import subprocess
import sys
import time
import webbrowser
from pathlib import Path


PORT = 5000
URL = f"http://localhost:{PORT}"


def debug(message):
    print(message, flush=True)


def wait_for_port(host="127.0.0.1", port=5000, timeout=60):
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.create_connection((host, port), timeout=2):
                return True
        except OSError:
            time.sleep(1)
    return False


def tail_log(log_file, max_chars=4000):
    if not log_file.exists():
        return ""
    try:
        return log_file.read_text(encoding="utf-8", errors="replace")[-max_chars:]
    except OSError as exc:
        return f"Could not read Streamlit log: {exc}"


def main():
    try:
        if getattr(sys, "frozen", False):
            app_dir = Path(sys.executable).resolve().parent
        else:
            app_dir = Path(__file__).resolve().parent

        log_dir = app_dir / "logs"
        log_dir.mkdir(exist_ok=True)
        log_file = log_dir / "corromics_streamlit.log"

        os.environ["ENVIRONMENT_MODE"] = "local"
        if os.name == "nt":
            os.environ["CORROMICS_GEMELLI_BACKEND"] = "disabled"
            debug("joint-RPCA/Gemelli: disabled for Windows executable")
        else:
            os.environ["CORROMICS_GEMELLI_BACKEND"] = "local"
            debug("joint-RPCA/Gemelli: enabled for local Linux/WSL environment")

        if os.name == "nt":
            streamlit_exe = app_dir / "envs" / "main" / "Scripts" / "streamlit.exe"
        else:
            streamlit_exe = app_dir / "envs" / "main" / "bin" / "streamlit"

        candidates = [
            app_dir / "Home.py",
            app_dir / "app.py",
            app_dir / "main.py",
            app_dir / "streamlit_app.py",
        ]

        app_file = next((candidate for candidate in candidates if candidate.exists()), None)
        if app_file is None:
            raise FileNotFoundError(
                "Could not find Streamlit entry file. Expected one of: "
                "Home.py, app.py, main.py, streamlit_app.py"
            )

        cmd = [
            str(streamlit_exe),
            "run",
            str(app_file),
            "--server.port",
            str(PORT),
            "--server.address",
            "127.0.0.1",
            "--server.headless",
            "true",
        ]

        debug(f"App directory: {app_dir}")
        debug(f"Selected Streamlit executable: {streamlit_exe}")
        debug(f"Selected app file: {app_file}")
        debug(f"Log file path: {log_file}")
        debug(f"Streamlit launch command: {' '.join(cmd)}")

        with open(log_file, "w", encoding="utf-8") as log:
            process = subprocess.Popen(
                cmd,
                cwd=str(app_dir),
                stdout=log,
                stderr=log,
                shell=False,
            )

        if wait_for_port(port=PORT):
            debug(f"Opening browser: {URL}")
            webbrowser.open(URL)
        else:
            debug("Streamlit did not start within 60 seconds.")
            debug("Last 4000 characters of Streamlit log:")
            debug(tail_log(log_file))
            input("Press Enter to close...")
            return 1

        return process.wait()
    except Exception as exc:
        debug(f"Startup failed: {exc}")
        input("Press Enter to close...")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
