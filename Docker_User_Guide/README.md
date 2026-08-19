# Run Corromics with Docker Desktop

This guide is for users who want to run Corromics locally with Docker on macOS
or Windows. Docker keeps the application and its dependencies inside a container,
so users do not need to install Python, Streamlit, or Gemelli manually.

The Corromics Docker image includes the Joint-RPCA/Gemelli worker environment.

## Files used by local Docker users

The required files are already included in the root of the Corromics repository:

- `Dockerfile` builds the application image.
- `run_server.sh` starts Streamlit on port 5000 inside the container.
- `environment-gemelli-worker.yml` defines the Joint-RPCA/Gemelli environment
  installed during the build.
- `requirements.txt` contains the main application dependencies.

The following files are intended for maintainers of the hosted GNPS2 deployment
and are not needed for normal local Docker use:

- `docker-compose.yml`
- `docker-compose-dev.yml`

## Hardware recommendations

- At least 8 GB of system RAM for small analyses.
- 16 GB RAM is recommended for normal local use.
- 32 GB or more is recommended for large datasets, SparCC, or Joint-RPCA.
- Keep enough free disk space for Docker Desktop and the Corromics image.

Docker Desktop itself supports systems with less memory, but the Corromics
analysis workload may require substantially more. Close other memory-intensive
applications during large runs.

## macOS setup

### 1. Install Docker Desktop

Download Docker Desktop using the official
[Docker Desktop for Mac installation guide](https://docs.docker.com/desktop/setup/install/mac-install/).

Choose the Docker Desktop installer that matches the Mac chip.

Open the downloaded `.dmg`, drag Docker into `Applications`, and start Docker
Desktop. Complete the initial setup prompts and wait until Docker reports that
the engine is running.

Verify the installation in Terminal:

```bash
docker --version
docker info
```

### 2. Clone Corromics

```bash
git clone https://github.com/abzer005/Corromics.git
cd Corromics
```

### 3. Build the image

#### Apple silicon Mac

The current Corromics Dockerfile installs an `x86_64` Linux environment. Build
the image with Docker's AMD64 emulation:

```bash
docker build --platform linux/amd64 -t corromics .
```

#### Intel Mac

```bash
docker build -t corromics .
```

The first build downloads the base image and all dependencies, so it can take
some time.

### 4. Run Corromics

#### Apple silicon Mac

```bash
docker run --rm \
  --platform linux/amd64 \
  -p 5000:5000 \
  -e ENVIRONMENT_MODE=local \
  corromics
```

#### Intel Mac

```bash
docker run --rm \
  -p 5000:5000 \
  -e ENVIRONMENT_MODE=local \
  corromics
```

Open <http://localhost:5000>. Keep the Terminal window open while using the
application. Press `Ctrl+C` to stop it.

## Windows setup

### 1. Prepare WSL2

Docker Desktop uses the WSL2 backend for Linux containers. Open PowerShell as
Administrator and run:

```powershell
wsl --install
wsl --update
```

Restart Windows if requested.

### 2. Install Docker Desktop

Follow the official
[Docker Desktop for Windows installation guide](https://docs.docker.com/desktop/setup/install/windows-install/).

During or after installation:

1. Start Docker Desktop.
2. Open **Settings → General**.
3. Enable **Use the WSL 2 based engine** if it is not already enabled.
4. Ensure Docker is using Linux containers.
5. Wait until Docker reports that the engine is running.

More details are available in Docker's
[WSL2 backend guide](https://docs.docker.com/desktop/features/wsl/).

Verify the installation in PowerShell:

```powershell
docker --version
docker info
```

### 3. Clone Corromics

If Git is not installed, install it from PowerShell:

```powershell
winget install --id Git.Git -e --source winget
```

Reopen PowerShell after installation and confirm that Git is available:

```powershell
git --version
```

Then clone Corromics:

```powershell
git clone https://github.com/abzer005/Corromics.git
cd Corromics
```

### 4. Build the image

```powershell
docker build -t corromics .
```

### 5. Run Corromics

```powershell
docker run --rm -p 5000:5000 -e ENVIRONMENT_MODE=local corromics
```

Open <http://localhost:5000>. Keep PowerShell open while using Corromics. Press
`Ctrl+C` to stop it.

## Optional: keep application logs on the host

Most users can skip this section.

The commands above start a temporary container using `--rm`. When that container
stops, Docker deletes it automatically, including any Corromics log files saved
inside it. This does not delete analysis results that users already downloaded
from the application.

Keeping the log files can be useful when troubleshooting an error or sharing
diagnostic information with the Corromics developers. The `-v` option below
connects the `logs` folder in the downloaded Corromics repository to the
container's `/app/logs` folder. Log files then remain on the computer after the
container stops.

macOS Terminal:

```bash
docker run --rm \
  -p 5000:5000 \
  -e ENVIRONMENT_MODE=local \
  -v "$(pwd)/logs:/app/logs" \
  corromics
```

Windows PowerShell:

```powershell
docker run --rm -p 5000:5000 -e ENVIRONMENT_MODE=local -v "${PWD}/logs:/app/logs" corromics
```

Apple silicon users should also include `--platform linux/amd64`.

## Optional: run without keeping the terminal open

The earlier `docker run --rm ...` command keeps Corromics attached to the open
Terminal or PowerShell window. This is the simplest option for most users: keep
the window open while working and press `Ctrl+C` when finished.

Background mode is an alternative, not an additional required step. It starts
Corromics without occupying the terminal, so the terminal can be closed or used
for other commands while Corromics continues running. The container is given the
name `corromics` so it can be stopped easily later.

```bash
docker run -d \
  --name corromics \
  -p 5000:5000 \
  -e ENVIRONMENT_MODE=local \
  corromics
```

When finished using Corromics, stop the background container and remove it:

```bash
docker stop corromics
docker rm corromics
```

Windows PowerShell should use the equivalent one-line command:

```powershell
docker run -d --name corromics -p 5000:5000 -e ENVIRONMENT_MODE=local corromics
```

Apple silicon users should add `--platform linux/amd64` to the `docker run`
command.

## Update Corromics

For a Git clone:

```bash
git pull
docker build --no-cache -t corromics .
```

Apple silicon users should retain `--platform linux/amd64` when rebuilding.

## Troubleshooting

### Docker is not running

Start Docker Desktop and wait until the engine is ready. Then check:

```bash
docker info
```

### Port 5000 is already in use

Use a different host port, such as 8502:

```bash
docker run --rm -p 8502:5000 -e ENVIRONMENT_MODE=local corromics
```

Open <http://localhost:8502>.

### The build or analysis runs out of memory

- Close other memory-intensive applications.
- Increase the resources available to Docker Desktop where supported.
- Reduce the number of omics features using filtering or hierarchical binning.
- Filter features before large SparCC analyses.
- Restart Docker Desktop after changing its resources.

Docker's current resource controls differ between macOS and the Windows WSL2
backend. See the official
[Docker Desktop settings guide](https://docs.docker.com/desktop/settings-and-maintenance/settings/).

### Apple silicon reports a platform mismatch

Ensure both build and run commands include:

```text
--platform linux/amd64
```

Running through emulation may be slower than running on an Intel/AMD64 machine.

### Windows Docker or WSL2 does not start

In PowerShell as Administrator:

```powershell
wsl --status
wsl --update
wsl --shutdown
```

Restart Docker Desktop afterward. Confirm **Use the WSL 2 based engine** is
enabled and Docker is using Linux containers.

### View container logs

For a background container:

```bash
docker logs corromics
```

Follow new log messages continuously:

```bash
docker logs -f corromics
```

## Docker Desktop licensing

Docker Desktop has its own license terms. Some larger organizations require a
paid Docker subscription for commercial use. Review the current terms in the
official Docker Desktop installation documentation before organizational use.
