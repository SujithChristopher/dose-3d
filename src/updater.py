"""Background check against GitHub Releases with an opt-in download-and-swap update.

Only active in the frozen (PyInstaller) exe - a dev checkout has no exe asset
to install, so the check is skipped when running from source. The network
check runs on a QThread so a slow or unreachable network never blocks the UI;
failures are swallowed silently rather than shown to the user.
"""
import json
import os
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Qt, Signal
from PySide6.QtWidgets import QMessageBox, QProgressDialog

from src._version import __version__

GITHUB_REPO = "SujithChristopher/dose-3d"
API_URL = f"https://api.github.com/repos/{GITHUB_REPO}/releases/latest"
REQUEST_TIMEOUT = 4  # seconds - don't block startup on a slow/offline network
ASSET_NAME = "EPID_Reconstruction_GUI.exe"


def _parse_version(v: str) -> tuple[int, ...]:
    v = v.lstrip("vV")
    parts = []
    for p in v.split("."):
        digits = "".join(c for c in p if c.isdigit())
        parts.append(int(digits) if digits else 0)
    return tuple(parts)


class _UpdateCheckWorker(QObject):
    found = Signal(str, str)  # version, download_url
    not_found = Signal()  # checked OK, already on the latest version
    check_error = Signal(str)  # message
    finished = Signal()

    def run(self):
        try:
            req = urllib.request.Request(
                API_URL, headers={"Accept": "application/vnd.github+json"}
            )
            with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as resp:
                data = json.load(resp)

            latest_tag = data.get("tag_name", "")
            if _parse_version(latest_tag) > _parse_version(__version__):
                download_url = next(
                    (a.get("browser_download_url") for a in data.get("assets", [])
                     if a.get("name") == ASSET_NAME),
                    None,
                )
                if download_url:
                    self.found.emit(latest_tag, download_url)
                else:
                    self.check_error.emit(f"Release {latest_tag} has no {ASSET_NAME} asset.")
            else:
                self.not_found.emit()
        except (urllib.error.URLError, ValueError, OSError) as e:
            self.check_error.emit(str(e))
        finally:
            self.finished.emit()


def check_for_update_async(parent):
    """Kick off a silent background update check on startup; prompts the user
    only if a newer release with a matching exe asset is found. Any other
    outcome (no update, network error) is intentionally not shown - use
    check_for_update_manual for a check that always reports its result."""
    if not getattr(sys, "frozen", False):
        return

    thread = QThread(parent)
    worker = _UpdateCheckWorker()
    worker.moveToThread(thread)
    thread.started.connect(worker.run)
    worker.found.connect(lambda version, url: _prompt_update(parent, version, url))
    worker.finished.connect(thread.quit)
    thread.finished.connect(thread.deleteLater)
    parent._update_thread = thread  # keep a reference so it isn't GC'd mid-flight
    thread.start()


def check_for_update_manual(parent):
    """User-triggered update check (e.g. Help menu) - always shows a dialog,
    whether that's an available update, "already up to date", or an error."""
    thread = QThread(parent)
    worker = _UpdateCheckWorker()
    worker.moveToThread(thread)
    thread.started.connect(worker.run)
    worker.found.connect(lambda version, url: _prompt_update(parent, version, url))
    worker.not_found.connect(lambda: QMessageBox.information(
        parent, "No updates available", f"You're running the latest version ({__version__})."))
    worker.check_error.connect(lambda msg: QMessageBox.warning(
        parent, "Update check failed", f"Could not check for updates:\n{msg}"))
    worker.finished.connect(thread.quit)
    thread.finished.connect(thread.deleteLater)
    parent._update_thread = thread  # keep a reference so it isn't GC'd mid-flight
    thread.start()


def _prompt_update(parent, version, download_url):
    if not getattr(sys, "frozen", False):
        QMessageBox.information(
            parent, "Update available",
            f"A new version ({version}) is available. You're running {__version__} from source.\n\n"
            f"Download it from: https://github.com/{GITHUB_REPO}/releases/latest")
        return

    reply = QMessageBox.question(
        parent,
        "Update available",
        f"A new version ({version}) is available. You're running {__version__}.\n\n"
        "Download and install it now? The app will restart.",
        QMessageBox.Yes | QMessageBox.No,
        QMessageBox.No,
    )
    if reply == QMessageBox.Yes:
        _download_and_swap(parent, download_url)


def _download_and_swap(parent, download_url):
    current_exe = Path(sys.executable)
    new_exe = Path(tempfile.gettempdir()) / f"{current_exe.stem}_new.exe"

    progress = QProgressDialog("Downloading update...", "Cancel", 0, 100, parent)
    progress.setWindowModality(Qt.WindowModal)
    progress.setMinimumDuration(0)
    progress.setValue(0)

    def _report(block_num, block_size, total_size):
        if total_size > 0:
            progress.setValue(min(100, int(block_num * block_size * 100 / total_size)))

    try:
        urllib.request.urlretrieve(download_url, new_exe, reporthook=_report)
    except (urllib.error.URLError, OSError) as e:
        progress.close()
        QMessageBox.warning(parent, "Update failed", f"Could not download update:\n{e}")
        return

    progress.close()

    helper_script = _write_helper_script(current_exe, new_exe)
    subprocess.Popen(
        ["powershell", "-NoProfile", "-WindowStyle", "Hidden", "-ExecutionPolicy", "Bypass",
         "-File", str(helper_script), str(os.getpid())],
        creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP,
        close_fds=True,
    )
    parent.close()


def _write_helper_script(current_exe: Path, new_exe: Path) -> Path:
    """A copy of the exe can't overwrite itself while running, so a detached
    helper waits for this process to exit, swaps the file, and relaunches."""
    script_path = Path(tempfile.gettempdir()) / "epid_gui_update.ps1"
    script = f"""
param([int]$ParentPid)
Wait-Process -Id $ParentPid -ErrorAction SilentlyContinue
Start-Sleep -Seconds 1
Copy-Item -Path "{new_exe}" -Destination "{current_exe}" -Force
Remove-Item -Path "{new_exe}" -Force -ErrorAction SilentlyContinue
Start-Process -FilePath "{current_exe}"
Remove-Item -Path $MyInvocation.MyCommand.Path -Force -ErrorAction SilentlyContinue
"""
    script_path.write_text(script, encoding="utf-8")
    return script_path
