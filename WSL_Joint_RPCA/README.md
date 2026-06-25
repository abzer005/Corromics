# Corromics joint-RPCA/Gemelli with WSL/Linux

Having WSL installed on Windows is useful beyond Corromics. Many bioinformatics and computational biology tools are developed primarily for Linux, or have dependencies that are more reliable on Linux than on native Windows. WSL lets Windows users run these Linux-compatible tools without needing a separate Linux computer.

The normal native Windows installation does not automatically use the WSL/Linux environment. The Windows installation remains available for the rest of the Corromics app, with joint-RPCA disabled.

To use joint-RPCA, run Corromics directly inside WSL/Linux with Gemelli and its required dependencies installed there.

This guide is written for users who are new to the terminal. In the examples below, the fake Ubuntu username is `maya`. Your username will be different.

## Section 1: If you already have WSL/Ubuntu installed

1. Open Ubuntu from the Start menu or Windows Terminal.

You should see a terminal window. The prompt may look something like this:

```text
maya@DESKTOP-12345:~$
```

This means you are inside Ubuntu. In this example, the Ubuntu username is `maya`.

2. Confirm Ubuntu/WSL is working:

```bash
uname -a
```

This should print a line of system information. It is okay if it looks technical.

3. Update Ubuntu:

```bash
sudo apt update
sudo apt upgrade -y
```

Ubuntu may ask for your password. Type the Ubuntu password you created when setting up WSL, then press Enter. The password will not visibly appear while you type. That is normal.

4. Clone or copy the Corromics repository into your Ubuntu home folder.

First go to your Ubuntu home folder:

```bash
cd ~
```

For a user named `maya`, this folder is:

```text
/home/maya
```

Then download the Corromics repository:

```bash
git clone https://github.com/Functional-Metabolomics-Lab/Corromics.git
cd Corromics
```

Use the Ubuntu home folder, such as `/home/maya/Corromics`, instead of putting the repository in your Windows Desktop or Downloads folder.

Why this matters: Ubuntu can also see Windows files, but those files appear under paths such as `/mnt/c/Users/maya/Desktop`. Beginners do not need to use those paths here. Keeping Corromics inside the Ubuntu home folder is usually faster and avoids Windows/Linux file permission problems.

5. Run the WSL/Linux joint-RPCA setup script.

The script may ask for your Ubuntu password if it needs to install missing Ubuntu packages such as `wget`. This is the Ubuntu password you created when setting up WSL, not your Windows password.

```bash
bash WSL_Joint_RPCA/setup_wsl_joint_rpca.sh
```

6. Start Corromics from WSL/Linux:

```bash
conda activate corromics-main
streamlit run Home.py --server.port 5000 --server.address 127.0.0.1
```

Keep this Ubuntu terminal window open while using Corromics. If you close it, the app will stop.

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

The native Windows installation and the WSL/Linux setup are separate paths.

Launching Corromics from the Windows installer shortcut starts the normal native Windows app. In that mode, joint-RPCA/Gemelli is disabled by design.

For joint-RPCA, start Corromics from the Ubuntu/WSL terminal after running `WSL_Joint_RPCA/setup_wsl_joint_rpca.sh`.
