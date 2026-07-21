# IT & Lab Administration Guide: Hosting and Exposing the Isoform Dashboard

This guide provides network administrators, IT managers, and lab managers with the system requirements, architectural patterns, and network configurations necessary to host the Isoform Dashboard locally and expose it securely to lab members.

---

## 1. Architectural & Deployment Overview

The Isoform Dashboard is a Python-based Dash application for visualizing genomic isoform distribution and 3D structural alignments (AlphaFold). Because genomic datasets and 3D structure folders are typically multi-gigabyte (or even terabyte) in size, the dashboard supports a **Compute/Storage Separation** architecture:
*   **Host Machine (Compute Server):** Runs the Python Dash backend, rendering interactive scatter plots and 3D protein structures using CPU/GPU resources.
*   **Client Machine (Researcher's Workstation):** Hosts the raw, large genomic datasets and accesses the dashboard via a web browser.
*   **Network Mapping:** Rather than transferring large files to the server or running resource-intensive visualization tools on local laptops, the client's data folder is mounted/shared over the local network (via SMB, CIFS, or SSHFS) so the host can read it directly.

```
┌────────────────────────────────────────┐               ┌────────────────────────────────────────┐
│     Client Workstation (Researcher)    │               │         Host Machine (Server)          │
│   (Houses large data, web browser)     │               │          (Runs Python Backend)         │
├────────────────────────────────────────┤               ├────────────────────────────────────────┤
│ [Data Folder]                          │ ◄───────────► │ [Network Mount Point]                  │
│  - GTF files, Fasta, 3D structures     │   SMB/CIFS    │  - Mounted as local drive or directory │
│                                        │   or SSHFS    │  - App reads files dynamically         │
│ [Web Browser]                          │ ◄───────────► │ [Dash Server Port: 8050]               │
│  - Accesses dashboard GUI              │  SSH Tunnel   │  - Generates UI layout & 3D models     │
│                                        │  or LAN IP    │                                        │
└────────────────────────────────────────┘               └────────────────────────────────────────┘
```

---

## 2. IT & Infrastructure Requirements

Before deploying the dashboard, ensure the network and server meet these specifications:

### System Requirements
*   **Operating System:** Windows 10/11, Windows Server, or Unix-like environments (Linux, macOS).
*   **Python:** Python 3.8 to 3.11.
*   **Memory (RAM):** 
    *   *Minimum:* 8 GB RAM (sufficient for small genomes or single-isoform datasets).
    *   *Recommended:* 16 GB+ RAM (required when loading large mammalian GTF files and complex AlphaFold geometry structures in-memory).

### Network & Security Requirements
*   **Port Requirements:**
    *   **Port `8050` (TCP):** Default port for the Dash web interface.
    *   **Port `22` (TCP):** Required if choosing secure SSH tunneling/mounting.
    *   **Port `445` (TCP):** Required for SMB/CIFS sharing (Windows-to-Windows or Windows-to-Linux).
*   **Access Control:** By default, the application has two network binding options:
    *   **Loopback Interface (`127.0.0.1`):** Restricts access to localhost. Highly secure. Combined with SSH port forwarding, it restricts dashboard access to authenticated SSH users only.
    *   **All Interfaces (`0.0.0.0`):** Exposes the dashboard to any machine on the Local Area Network (LAN). Requires network firewall clearance.

---

## 3. Host Server Environment Preparation

Follow these steps to prepare the Python environment on the host machine.

### Option A: Windows Host/Server Setup
1. **Clone the Repository:**
   ```powershell
   git clone https://github.com/borbalamuszka/isoform-analysis.git
   cd isoform-analysis
   ```
2. **Configure Virtual Environment:**
   Initialize and activate a virtual environment to isolate python dependencies:
   ```powershell
   python -m venv isoform_dashboard_env
   isoform_dashboard_env\Scripts\activate
   ```
3. **Install Dependencies:**
   ```powershell
   pip install -r requirements.txt
   ```

### Option B: Linux/Unix Host/Server Setup
1. **Clone the Repository:**
   ```bash
   git clone https://github.com/borbalamuszka/isoform-analysis.git
   cd isoform-analysis
   ```
2. **Configure Virtual Environment:**
   ```bash
   python3 -m venv isoform_dashboard_env
   source isoform_dashboard_env/bin/activate
   ```
3. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

---

## 4. Storage Integration & Data Mounting

To access data stored on client workstations, you must map their directories to the server. Choose the protocol that aligns with your lab's operating system environment.

### Protocol A: SMB (Windows Client to Windows Host)
*Best for pure Windows environments. Utilizes native Windows sharing.*

1. **On the Client Workstation (Storage Origin):**
   *   Right-click the data folder (e.g., `C:\IsoformData`) and select **Properties**.
   *   Go to **Sharing** -> **Advanced Sharing...** and check **Share this folder**.
   *   Click **Permissions**. Set the permissions to **Read** for the network user connecting from the server.
   *   Note the client's IP address (e.g., `192.168.1.15`) and the share name (e.g., `IsoformData`).
2. **On the Host Server (Mount Target):**
   *   Open PowerShell and mount the shared directory as a virtual network drive (e.g., `Z:`):
     ```powershell
     net use Z: \\192.168.1.15\IsoformData /persistent:yes
     ```
     *(Provide client-side credentials when prompted).*

### Protocol B: CIFS (Windows Client to Linux Host)
*Allows a Linux compute server to read data hosted on a Windows workstation.*

1. **On the Client Workstation (Windows):**
   *   Configure the SMB share as described in Protocol A.
2. **On the Host Server (Linux):**
   *   Install CIFS utilities:
     ```bash
     sudo apt-get update && sudo apt-get install cifs-utils -y
     ```
   *   Create a local mount point:
     ```bash
     sudo mkdir -p /mnt/client_data
     ```
   *   Mount the share:
     ```bash
     sudo mount -t cifs -o username=<Windows_Username>,password=<Windows_Password>,iocharset=utf8,file_mode=0777,dir_mode=0777 //192.168.1.15/IsoformData /mnt/client_data
     ```

### Protocol C: SSHFS (Unix/macOS Client to Unix Host)
*Secure file-system mounting over SSH.*

1. **On the Host Server (Linux/macOS):**
   *   Install SSHFS:
       *   **Debian/Ubuntu:** `sudo apt install sshfs`
       *   **macOS:** `brew install macfuse && brew install gromgit/fuse/sshfs`
   *   Mount the remote client directory (assumes client workstation has SSH service running):
     ```bash
     mkdir -p ~/client_data
     sshfs client_user@192.168.1.15:/path/to/client/data ~/client_data
     ```

---

## 5. Network Access & Exposure Options

Lab managers and IT administrators can expose the dashboard using two models depending on security requirements.

### Option 1: Secure SSH Tunneling (Recommended for Confidential Data)
Exposing the application on `127.0.0.1` ensures that only users with valid SSH accounts on the host machine can access the dashboard. This prevents unauthorized users on the local network from viewing proprietary genomic designs.

```
[Browser] ──► (Port 8050) ──► [SSH Client] ── Encrypted Tunnel ──► [SSH Daemon] ──► (Port 8050) ──► [Dash App]
   ├─────────────────────────────────────┤                           ├──────────────────────────────────┤
   │           Client Machine            │                           │           Host Server            │
```

1.  **On the Host Server:** Run the dashboard bound to the loopback interface:
    ```bash
    python -m isoform_dashboard.dashboard_app --host 127.0.0.1 --port 8050 [arguments...]
    ```
2.  **On the Client Workstation:** Forward traffic on the local port `8050` through SSH to the server:
    ```bash
    ssh -L 8050:127.0.0.1:8050 host_user@192.168.1.100
    ```
3.  **Client Browsing:** The researcher accesses the interface locally at `http://127.0.0.1:8050`.

---

### Option 2: Local Network Binding (Direct Access)
Exposing the server on `0.0.0.0` allows direct browser access via the host server's LAN IP address. This is convenient for open environments but requires open firewall ports.

1.  **On the Host Server:** Run the dashboard bound to all interfaces:
    ```bash
    python -m isoform_dashboard.dashboard_app --host 0.0.0.0 --port 8050 [arguments...]
    ```
2.  **Firewall Administration:**
    *   **Windows Host Firewall:** Run in administrative PowerShell:
        ```powershell
        New-NetFirewallRule -DisplayName "Isoform Dashboard Port 8050" -Direction Inbound -Action Allow -Protocol TCP -LocalPort 8050
        ```
    *   **Linux Host Firewall (`ufw`):**
        ```bash
        sudo ufw allow 8050/tcp
        ```
3.  **Client Browsing:** Client workstations access the dashboard using the host server's IP address:
    `http://<host-ip-address>:8050` (e.g., `http://192.168.1.100:8050`).

---

## 6. End-to-End Walkthrough (Windows-to-Windows deployment)

Here is a deployment template checklist for a common Windows-based lab setup:

- [ ] **Step 1 (Client Workstation - IP `192.168.1.15`):**
  Share the folder `C:\IsoformData` under the share name `IsoformData` (Read permissions for `LabServerUser`).
- [ ] **Step 2 (Host Server - IP `192.168.1.100`):**
  Mount the share to the `Z:` drive:
  ```powershell
  net use Z: \\192.168.1.15\IsoformData /persistent:yes
  ```
- [ ] **Step 3 (Host Server):**
  Open PowerShell, navigate to the dashboard directory, activate the environment, and launch the service:
  ```powershell
  cd C:\Users\vilsn\Documents\isoform-analysis
  isoform_dashboard_env\Scripts\activate
  python -m isoform_dashboard.dashboard_app `
    --host 127.0.0.1 `
    --port 8050 `
    --input-mean Z:\distributions_mean.tsv `
    --input-sum Z:\distributions_sum.tsv `
    --ci-file Z:\confidence_intervals.tsv `
    --exons Z:\expressed_isoforms.gtf `
    --proteins Z:\proteins.fasta `
    --geometry-dir Z:\alphafold_geometry
  ```
- [ ] **Step 4 (Client Workstation):**
  Establish the secure port forwarding session:
  ```powershell
  ssh -L 8050:127.0.0.1:8050 LabServerUser@192.168.1.100
  ```
  *(Keep this terminal open while using the dashboard).*
- [ ] **Step 5 (Client Workstation):**
  Navigate to `http://127.0.0.1:8050` in a web browser.

---

## 7. Dynamic Data Reloading & Administrative Considerations

To simplify maintenance, the dashboard supports a **zero-configuration launch state**. This allows you to run the dashboard as a persistent background service without hardcoding startup file paths.

### Zero-Configuration Service Launch
Start the service with no input arguments specified:
```powershell
python -m isoform_dashboard.dashboard_app --host 0.0.0.0 --port 8050
```
Upon connection, users are greeted with a configuration banner.

### Graphical Path Selector & Server-Side Security
Through the **⚙ Configure Data Sources** dashboard utility, users can select files and folders dynamically:
1.  **Browser Interface:** Users click **Browse...** next to any input field.
2.  **Server File Navigation:** A file selection window displays files and directories *on the host filesystem*.
3.  **Data Mounting Detection:** Mapped network drives (such as `Z:\`) or mount directories (like `/mnt`) are automatically detected and available under "Quick Jump" buttons.
4.  **Automatic In-Memory Recomputation:** Once files are selected, clicking **Apply Changes & Reload** reads the files, performs calculations in-memory, and updates the visualization state. **No server restarts are required.**

> [!WARNING]
> **IT Security Note:** The in-app visual file browser navigates the *host machine's filesystem*. To prevent unauthorized directory traversal, ensure the user account running the Python Dash application has restricted permissions on the host system. It should only have read access to target data folders and mount points, and no read/write access to sensitive system files.
