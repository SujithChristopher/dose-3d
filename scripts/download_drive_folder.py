"""Download a public Google Drive folder (recursively) into a local directory.

Used to fetch EPID DICOM datasets that are shared as Drive folders instead of
being checked into the repo.

Usage:
    uv run python scripts/download_drive_folder.py <folder-url-or-id> -o dataset/my-case
    uv run python scripts/download_drive_folder.py <folder-url-or-id> --dry-run

The folder must be shared as "Anyone with the link". Existing files are skipped,
so an interrupted run can simply be re-run to resume.

Both phases are parallel and time-bounded:

* Listing walks the folder tree with a thread pool, one request per folder, and
  caches the result to a manifest JSON next to the output so re-runs skip it.
  (gdown's own walker is serial and issues requests with no timeout, so a single
  stalled connection hangs the whole listing.)
* Downloading streams files over keep-alive sessions, one per worker thread --
  no repeated TLS handshake, no HTML round-trip for files small enough to skip
  Drive's virus-scan warning. gdown is the fallback for anything it refuses.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import threading
import time
import urllib.parse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
from pathlib import Path

try:
    import bs4
    import gdown
    import requests
    from gdown.download import _sanitize_filename
    from requests.adapters import HTTPAdapter
except ImportError:  # pragma: no cover - dependency hint
    sys.exit("gdown is required: uv add gdown  (or: uv run --with gdown python ...)")

from tqdm import tqdm

DEFAULT_URL = "https://drive.google.com/drive/folders/17Y373l-n7C6IZy7JAPZVxJHT05Q7hy7f?usp=drive_link"
DIRECT_URL = "https://drive.usercontent.google.com/download"
FOLDER_VIEW_URL = "https://drive.google.com/embeddedfolderview"
MANIFEST_NAME = ".drive_manifest.json"
CHUNK = 1 << 20  # 1 MiB
LIST_TIMEOUT = (10, 30)  # (connect, read) seconds

FILE_HREF = re.compile(r"https://drive\.google\.com/file/d/([-\w]{25,})/view")
DOCS_HREF = re.compile(r"https://docs\.google\.com/\w+/d/([-\w]{25,})/")
FOLDER_HREF = re.compile(r"https://drive\.google\.com/drive/folders/([-\w]{25,})")

_local = threading.local()


def get_session(pool_size: int) -> requests.Session:
    """One keep-alive session per worker thread."""
    session = getattr(_local, "session", None)
    if session is None:
        session = requests.Session()
        adapter = HTTPAdapter(pool_connections=pool_size, pool_maxsize=pool_size)
        session.mount("https://", adapter)
        _local.session = session
    return session


def extract_folder_id(url: str) -> str:
    if "/" not in url:
        return url
    return urllib.parse.urlparse(url).path.rstrip("/").split("/")[-1]


def parse_folder(folder_id: str, retries: int, pool_size: int) -> list[tuple[str, str, bool]]:
    """Return one folder's children as (id, name, is_folder), non-recursively."""
    session = get_session(pool_size)
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            response = session.get(
                FOLDER_VIEW_URL, params={"id": folder_id}, timeout=LIST_TIMEOUT
            )
            response.raise_for_status()
            break
        except Exception as error:  # transient network / rate limit
            last_error = error
            if attempt + 1 < retries:
                time.sleep(2**attempt)
    else:
        raise RuntimeError(f"Could not list folder {folder_id}: {last_error}")

    soup = bs4.BeautifulSoup(response.text, features="html.parser")
    children: list[tuple[str, str, bool]] = []
    for a_tag in soup.find_all(name="a"):
        href = a_tag.get("href", "")
        if not isinstance(href, str):
            continue
        name = _sanitize_filename(filename=a_tag.get_text(strip=True))
        folder_match = FOLDER_HREF.match(href)
        if folder_match:
            children.append((folder_match.group(1), name, True))
            continue
        file_match = FILE_HREF.match(href) or DOCS_HREF.match(href)
        if file_match:
            children.append((file_match.group(1), name, False))
    return children


def list_folder(url: str, workers: int, retries: int) -> list[tuple[str, str]]:
    """Walk the folder tree in parallel. Returns (file_id, relative_path) pairs."""
    root_id = extract_folder_id(url)
    files: list[tuple[str, str]] = []
    seen_folders = {root_id}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(parse_folder, root_id, retries, workers): Path(".")}
        with tqdm(desc="Listing", unit="folder") as bar:
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    prefix = pending.pop(future)
                    for child_id, name, is_folder in future.result():
                        if is_folder:
                            if child_id in seen_folders:  # shortcut / cycle guard
                                continue
                            seen_folders.add(child_id)
                            pending[
                                pool.submit(parse_folder, child_id, retries, workers)
                            ] = prefix / name
                        else:
                            files.append((child_id, str(prefix / name)))
                    bar.update(1)
                    bar.set_postfix_str(f"{len(files)} files")
    if not files:
        raise RuntimeError(
            "No files found. Check that the link is valid and shared with "
            "'Anyone with the link'."
        )
    files.sort(key=lambda item: item[1])
    return files


def load_manifest(path: Path, url: str) -> list[tuple[str, str]] | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if data.get("url") != url:
        return None
    return [(entry["id"], entry["path"]) for entry in data.get("files", [])]


def save_manifest(path: Path, url: str, files: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "url": url,
        "listed_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "files": [{"id": file_id, "path": rel} for file_id, rel in files],
    }
    path.write_text(json.dumps(payload, indent=1))


def direct_download(file_id: str, part: Path, pool_size: int) -> int:
    """Stream one file straight from Drive. Returns bytes written, 0 if refused."""
    session = get_session(pool_size)
    params = {"id": file_id, "export": "download", "confirm": "t"}
    with session.get(DIRECT_URL, params=params, stream=True, timeout=(10, 60)) as response:
        response.raise_for_status()
        # A confirm/quota/permission page comes back as HTML instead of the file.
        if "text/html" in response.headers.get("Content-Type", ""):
            return 0
        written = 0
        with open(part, "wb") as handle:
            for chunk in response.iter_content(CHUNK):
                handle.write(chunk)
                written += len(chunk)
    return written


def fetch_one(file_id: str, dest: Path, retries: int, pool_size: int) -> tuple[Path, int, bool]:
    """Download a single Drive file, retrying transient failures.

    Returns (path, bytes, ok). Writes to a .part file and renames on success so
    an interrupted run never leaves a truncated file that a resume would skip.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            size = direct_download(file_id, part, pool_size)
        except Exception:
            size = 0
        if size:
            part.replace(dest)
            return dest, size, True
        part.unlink(missing_ok=True)
        # Direct URL refused (confirm page, transient error): let gdown handle it.
        try:
            if gdown.download(id=file_id, output=str(part), quiet=True) is not None:
                if part.exists() and part.stat().st_size > 0:
                    size = part.stat().st_size
                    part.replace(dest)
                    return dest, size, True
        except Exception:
            pass
        part.unlink(missing_ok=True)
        if attempt + 1 < retries:
            time.sleep(2**attempt)
    return dest, 0, False


def run(args: argparse.Namespace) -> int:
    manifest_path = args.manifest or args.output / MANIFEST_NAME
    files = None if args.refresh else load_manifest(manifest_path, args.url)
    if files is None:
        print(f"Listing {args.url} ...")
        files = list_folder(args.url, args.list_workers, args.retries)
        save_manifest(manifest_path, args.url, files)
        print(f"Found {len(files)} files (cached in {manifest_path})")
    else:
        print(f"Using cached listing from {manifest_path} ({len(files)} files, --refresh to re-list)")

    targets = [(file_id, args.output / rel) for file_id, rel in files]
    if args.limit:
        targets = targets[: args.limit]

    if args.dry_run:
        for _, path in targets:
            print(path)
        return 0

    pending = [
        (file_id, path)
        for file_id, path in targets
        if args.overwrite or not (path.exists() and path.stat().st_size > 0)
    ]
    skipped = len(targets) - len(pending)
    if skipped:
        print(f"Skipping {skipped} files already present")
    if not pending:
        print("Nothing to do.")
        return 0

    print(f"Downloading {len(pending)} files with {args.workers} workers ...")
    failed: list[Path] = []
    total_bytes = 0
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(fetch_one, file_id, path, args.retries, args.workers)
            for file_id, path in pending
        ]
        try:
            with tqdm(total=len(futures), unit="file") as bar:
                for future in as_completed(futures):
                    path, size, ok = future.result()
                    if ok:
                        total_bytes += size
                    else:
                        failed.append(path)
                    bar.update(1)
                    elapsed = time.monotonic() - started
                    bar.set_postfix_str(f"{total_bytes / 1e6 / max(elapsed, 1e-9):.1f} MB/s")
        except KeyboardInterrupt:
            pool.shutdown(wait=False, cancel_futures=True)
            raise

    elapsed = time.monotonic() - started
    print(
        f"Downloaded {len(pending) - len(failed)}/{len(pending)} files "
        f"({total_bytes / 1e6:.0f} MB in {elapsed:.0f}s, "
        f"{total_bytes / 1e6 / max(elapsed, 1e-9):.1f} MB/s) into {args.output}"
    )
    if failed:
        print(f"{len(failed)} failed (re-run to retry):")
        for path in failed[:20]:
            print(f"  {path}")
        if len(failed) > 20:
            print(f"  ... and {len(failed) - 20} more")
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("url", nargs="?", default=DEFAULT_URL, help="Drive folder URL or ID")
    parser.add_argument("-o", "--output", default=Path("dataset"), type=Path, help="local destination directory")
    parser.add_argument("--dry-run", action="store_true", help="list files only, download nothing")
    parser.add_argument("--workers", type=int, default=16, help="parallel downloads (default: 16)")
    parser.add_argument("--list-workers", type=int, default=8, help="parallel folder listings (default: 8)")
    parser.add_argument("--retries", type=int, default=3, help="attempts per request (default: 3)")
    parser.add_argument("--overwrite", action="store_true", help="re-download files that already exist")
    parser.add_argument("--limit", type=int, help="stop after N files (useful for a quick test)")
    parser.add_argument("--manifest", type=Path, help=f"listing cache path (default: <output>/{MANIFEST_NAME})")
    parser.add_argument("--refresh", action="store_true", help="ignore the cached listing and re-walk the folder")
    args = parser.parse_args()

    try:
        return run(args)
    except KeyboardInterrupt:
        print("\nInterrupted. Re-run to resume; completed files are kept.", file=sys.stderr)
        return 130
    except RuntimeError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
