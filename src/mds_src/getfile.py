#!/usr/bin/env python3
"""
getfile.py

Download a file given a Google Drive or OneDrive share URL.

Usage:
    python getfile.py <share_url> [dest_path]

If dest_path is omitted the filename 'OUTCAR_collect' will be used.
"""

from __future__ import annotations
import re
import requests
import sys
from urllib.parse import urlparse, parse_qs, urljoin

# ---------------- Google Drive helpers ----------------
def _extract_drive_file_id(url: str) -> str:
    """Extract Google Drive file id from common share URL patterns."""
    m = re.search(r"/file/d/([a-zA-Z0-9_-]+)", url)
    if m:
        return m.group(1)
    qs = parse_qs(urlparse(url).query)
    if "id" in qs:
        return qs["id"][0]
    # also support sharing urls that embed id near 'usp=sharing' style
    m2 = re.search(r"[-\w]{25,}", url)   # fallback: long token
    if m2:
        return m2.group(0)
    raise ValueError("Could not extract Google Drive file id")

def download_google_drive_file(share_url: str, dest_path: str, chunk_size: int = 1024 * 1024):
    """Download file from Google Drive public share link, handling the 'confirm' step."""
    file_id = _extract_drive_file_id(share_url)

    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    base_url = "https://drive.google.com/uc"
    params = {"export": "download", "id": file_id}
    res = session.get(base_url, params=params, stream=True)

    # If the response is HTML (e.g. a confirm page for large files), try to find a confirm token or redirect endpoint
    content_type = res.headers.get("Content-Type", "")
    if "text/html" in content_type:
        # try to extract hidden inputs from the HTML
        hidden_inputs = re.findall(r'<input[^>]+name="([^"]+)"\s+value="([^"]*)"', res.text)
        data = {k: v for k, v in hidden_inputs}
        if "confirm" in data:
            # Some Drive pages forward to a different endpoint; try the usercontent endpoint used in some flows
            confirm_url = "https://drive.usercontent.google.com/download"
            res = session.get(confirm_url, params=data, stream=True)
        else:
            # Another fallback: try the "download" URL that sometimes appears as a link
            m = re.search(r'href="([^"]*confirm=[^"]+)"', res.text)
            if m:
                redirect_url = m.group(1)
                if redirect_url.startswith("/"):
                    redirect_url = urljoin("https://drive.google.com", redirect_url)
                res = session.get(redirect_url, stream=True)
            else:
                raise RuntimeError("Google Drive: couldn't find confirm token — check file permissions (set 'Anyone with link').")

    res.raise_for_status()
    print(f"Downloading to {dest_path}...")
    with open(dest_path, "wb") as f:
        for chunk in res.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    print("Download complete.")
    return dest_path

# ---------------- OneDrive helpers ----------------
def _is_onedrive_url(parsed) -> bool:
    host = (parsed.netloc or "").lower()
    return any(x in host for x in ("1drv.ms", "onedrive.live.com", "sharepoint.com", "onedrive.com"))

def download_onedrive_file(share_url: str, dest_path: str, chunk_size: int = 1024 * 1024):
    """
    Download file from OneDrive/1drv.ms/sharepoint public share links.

    Strategy:
      1) Follow redirects (requests will follow short 1drv.ms links).
      2) If final response has Content-Disposition -> stream it.
      3) Otherwise try appending '?download=1' to final URL and request again.
      4) If still HTML, try to parse out a direct download link in the HTML (common patterns).
      5) Stream to disk.
    """
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    # initial request (follow redirects to reach the real landing page / download link)
    r = session.get(share_url, allow_redirects=True, stream=True)
    final_url = r.url
    # If response is already binary (attachment), stream it
    if "content-disposition" in r.headers or not r.headers.get("Content-Type", "").startswith("text/html"):
        r.raise_for_status()
        print(f"Downloading (direct) to {dest_path}...")
        with open(dest_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
        return dest_path

    # Try adding download=1 param to force download
    if "download=1" not in final_url:
        trial = final_url + ("&" if "?" in final_url else "?") + "download=1"
        r2 = session.get(trial, allow_redirects=True, stream=True)
        if "content-disposition" in r2.headers or not r2.headers.get("Content-Type", "").startswith("text/html"):
            r2.raise_for_status()
            print(f"Downloading (download=1) to {dest_path}...")
            with open(dest_path, "wb") as f:
                for chunk in r2.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
            print("Download complete.")
            return dest_path
        # otherwise keep the HTML for parsing
        r = r2
        final_url = r.url

    # If we reached here we have an HTML landing page — try to parse for direct download links
    text = r.text

    # Try common anchor patterns that point to download URLs
    # Examples:
    #   <a href="https://onedrive.live.com/download?cid=...&resid=...">Download</a>
    #   window.location = "https://public.bn.files.1drv.com/..."
    m = re.search(r'href=["\']([^"\']*download[^"\']*)["\']', text, re.IGNORECASE)
    if not m:
        m = re.search(r'window\.location(?:\.href)?\s*=\s*["\']([^"\']+)["\']', text, re.IGNORECASE)
    if not m:
        # try to find large direct links pointing to files (cdn / public)
        m = re.search(r'https?://[^\s"\']+\.1drv\.com[^\s"\']*', text)
    if m:
        dl = m.group(1)
        if dl.startswith("/"):
            dl = urljoin(final_url, dl)
        r3 = session.get(dl, stream=True)
        r3.raise_for_status()
        print(f"Downloading (extracted link) to {dest_path}...")
        with open(dest_path, "wb") as f:
            for chunk in r3.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
        return dest_path

    # If all parsing attempts fail, raise informative error
    raise RuntimeError("Could not locate direct OneDrive download link. "
                       "Ensure the file is shared publicly ('Anyone with the link') "
                       "or use an alternative sharing method.")

# ---------------- General entrypoint ----------------
def download_shared_file(share_url: str, dest_path: str = "OUTCAR_collect"):
    """
    Autodetect service and download accordingly.
    Supports Google Drive and OneDrive (public share links).
    """
    parsed = urlparse(share_url)
    host = (parsed.netloc or "").lower()
    if "drive.google" in host or "docs.google" in host or "googleusercontent.com" in host or "drive.google.com" in share_url:
        return download_google_drive_file(share_url, dest_path)
    if _is_onedrive_url(parsed):
        return download_onedrive_file(share_url, dest_path)
    # fallback: try to GET directly (handles direct HTTP links)
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})
    r = session.get(share_url, stream=True, allow_redirects=True)
    # If we get a binary response stream it
    if "content-disposition" in r.headers or not r.headers.get("Content-Type", "").startswith("text/html"):
        r.raise_for_status()
        print(f"Downloading (generic) to {dest_path}...")
        with open(dest_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                if chunk:
                    f.write(chunk)
        print("Download complete.")
        return dest_path

    # Otherwise error
    raise RuntimeError("Unknown share link type or can't find direct download. "
                       "Supported: Google Drive, OneDrive, or direct file URLs.")

# ---------------- CLI ----------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python getfile.py <share_url> [dest_path]")
        sys.exit(1)
    url = sys.argv[1]
    dest = sys.argv[2] if len(sys.argv) >= 3 else "OUTCAR_collect"
    try:
        out = download_shared_file(url, dest)
        print("Saved to:", out)
    except Exception as e:
        print("ERROR:", e)
        sys.exit(2)
