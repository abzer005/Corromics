# Corromics joint-RPCA/Gemelli with WSL/Linux

Having WSL installed on Windows is useful beyond Corromics. Many bioinformatics and computational biology tools are developed primarily for Linux, or have dependencies that are more reliable on Linux than on native Windows. WSL lets Windows users run these Linux-compatible tools without needing a separate Linux computer.

The normal Windows `.exe` does not automatically use the WSL/Linux environment. The Windows executable remains available for the rest of the Corromics app, with joint-RPCA disabled.

To use joint-RPCA, run Corromics directly inside WSL/Linux with Gemelli and its required dependencies installed there.

## Section 1: If you already have WSL/Ubuntu installed

1. Open Ubuntu from the Start menu or Windows Terminal.

2. Confirm WSL is working:

```bash
uname -a
```

3. Update Ubuntu:

```bash
sudo apt update
sudo apt upgrade -y
```

4. Clone or copy the Corromics repository into the Ubuntu/Linux filesystem, preferably under your home directory:

```bash
cd ~
git clone <CORROMICS_REPOSITORY_URL>
cd Corromics
```

Using the Ubuntu home folder is preferred over `/mnt/c/...` because it is usually faster and avoids some Windows/Linux file permission issues.

5. Run the WSL/Linux joint-RPCA setup script.

The script may ask for your Ubuntu password if it needs to install missing Ubuntu packages such as `wget`. This is the Ubuntu password you created when setting up WSL, not your Windows password.

```bash
bash setup_wsl_joint_rpca.sh
```

6. Start Corromics from WSL/Linux:

```bash
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

7. Open Corromics in the Windows browser:

```text
http://localhost:5000
```

## Section 2: If you do not have WSL/Ubuntu installed

1. Open PowerShell or Windows Terminal as Administrator.

2. Install Ubuntu with WSL:

```powershell
wsl --install -d Ubuntu
```

3. Restart Windows if prompted.

4. Open Ubuntu from the Start menu.

5. Create the Ubuntu username/password when prompted.

6. Then follow Section 1.

## Important Notes

The Windows `.exe` and the WSL/Linux setup are separate paths.

Double-clicking `Corromics.exe` on Windows starts the normal Windows executable. In that mode, joint-RPCA/Gemelli is disabled by design.

For joint-RPCA, start Corromics from the Ubuntu/WSL terminal after running `setup_wsl_joint_rpca.sh`.
