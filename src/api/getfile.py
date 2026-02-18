#!/usr/bin/env python3
"""
getfile.py

Download a file from a share link (Google Drive or Dropbox) and save as OUTCAR_collect
(or other dest path). Designed to be called from other SLUSCHI scripts.

Usage:
    python getfile.py "<share-url>" [dest_path]

If dest_path not provided, defaults to "OUTCAR_collect".
"""

from __future__ import annotations
import re
import sys
import requests
from urllib.parse import urlparse, parse_qs, urlencode, urlunparse

# ---------------- Google Drive helpers (your existing approach) ----------------
def _extract_drive_file_id(url: str) -> str:
    """Extract Google Drive file id from common share URL formats."""
    m = re.search(r"/file/d/([a-zA-Z0-9_-]+)", url)
    if m:
        return m.group(1)
    qs = parse_qs(urlparse(url).query)
    if "id" in qs:
        return qs["id"][0]
    # fallback: maybe the url is the id itself
    if re.match(r'^[a-zA-Z0-9_-]{10,}$', url):
        return url
    raise ValueError("Could not extract Google Drive file id from URL")

def download_google_drive_file(share_url: str, dest_path: str, chunk_size: int = 1024 * 1024):
    """
    Download a Google Drive shared file (handles the confirmation HTML page for large files).
    """
    file_id = _extract_drive_file_id(share_url)
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    base_url = "https://drive.google.com/uc"
    params = {"export": "download", "id": file_id}
    res = session.get(base_url, params=params, stream=True)

    # If Google returns an HTML page (large file confirmation), try to extract hidden fields
    if "text/html" in res.headers.get("Content-Type", ""):
        hidden_inputs = re.findall(r'name="([^"]+)"\s+value="([^"]+)"', res.text)
        data = {k: v for k, v in hidden_inputs}
        if "confirm" in data:
            # alternative endpoint (works with certain Google pages)
            confirm_url = "https://drive.usercontent.google.com/download"
            res = session.get(confirm_url, params=data, stream=True)
        else:
            raise RuntimeError("Google Drive returned an HTML page. Ensure 'Anyone with link' is enabled.")

    res.raise_for_status()
    print(f"Downloading Google Drive file to {dest_path} ...")
    with open(dest_path, "wb") as f:
        for chunk in res.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    print("Google Drive download complete.")
    return dest_path

# ---------------- Dropbox helpers ----------------
def _is_dropbox_url(url: str) -> bool:
    net = urlparse(url).netloc.lower()
    return "dropbox.com" in net or "dropboxusercontent.com" in net

def _dropbox_direct_url(share_url: str) -> str:
    """
    Convert common Dropbox share URLs to a direct-download link (dl=1).
    Examples:
      - https://www.dropbox.com/s/<id>/<name>?dl=0 -> change to dl=1
      - https://dl.dropboxusercontent.com/...    -> already direct
    """
    parsed = urlparse(share_url)
    net = parsed.netloc.lower()

    if "dropboxusercontent.com" in net:
        # already a direct content URL (return as-is)
        return share_url

    if "dropbox.com" in net:
        qs = parse_qs(parsed.query)
        qs["dl"] = ["1"]
        new_query = urlencode(qs, doseq=True)
        new_parsed = parsed._replace(query=new_query)
        return urlunparse(new_parsed)

    raise ValueError("Not a Dropbox share URL")

def download_dropbox_file(share_url: str, dest_path: str, chunk_size: int = 1024 * 1024):
    """Download a Dropbox shared file by converting share link to direct-transfer (dl=1)."""
    direct = _dropbox_direct_url(share_url)
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})
    res = session.get(direct, stream=True, allow_redirects=True)
    if res.status_code == 403:
        raise RuntimeError("Access forbidden: check Dropbox sharing permissions (anyone with link).")
    res.raise_for_status()
    print(f"Downloading Dropbox file to {dest_path} ...")
    with open(dest_path, "wb") as f:
        for chunk in res.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    print("Dropbox download complete.")
    return dest_path

# ---------------- Dispatcher ----------------
def download_from_link(share_url: str, dest_path: str = "OUTCAR_collect"):
    """Auto-detect provider (Google Drive or Dropbox) and download to dest_path."""
    parsed = urlparse(share_url)
    net = parsed.netloc.lower()
    if "drive.google.com" in net or "docs.google.com" in net or "googleusercontent.com" in net:
        return download_google_drive_file(share_url, dest_path)
    if _is_dropbox_url(share_url):
        return download_dropbox_file(share_url, dest_path)
    raise ValueError("Unsupported provider. URL must be Google Drive or Dropbox share link.")

# ---------------- CLI ----------------
def main():
    if len(sys.argv) < 2:
        print("Usage: python getfile.py <share-url> [dest_path]")
        raise SystemExit(1)
    url = sys.argv[1].strip()
    dest = sys.argv[2].strip() if len(sys.argv) >= 3 else "OUTCAR_collect"
    try:
        print("Requested URL:", url)
        download_from_link(url, dest)
    except Exception as e:
        print("ERROR:", e)
        raise

if __name__ == "__main__":
    main()
