import re
import requests
from urllib.parse import urlparse, parse_qs

def _extract_drive_file_id(url: str) -> str:
    m = re.search(r"/file/d/([a-zA-Z0-9_-]+)", url)
    if m: return m.group(1)
    qs = parse_qs(urlparse(url).query)
    if "id" in qs: return qs["id"][0]
    raise ValueError("Could not extract Google Drive file id")

def download_google_drive_file(share_url: str, dest_path: str, chunk_size: int = 1024 * 1024):
    file_id = _extract_drive_file_id(share_url)
    
    # Use a session to automatically handle cookies
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0"})

    # Step 1: Initial Request
    base_url = "https://drive.google.com/uc"
    params = {"export": "download", "id": file_id}
    res = session.get(base_url, params=params, stream=True)

    # Step 2: Check for the "Virus Scan Warning" page
    if "text/html" in res.headers.get("Content-Type", ""):
        # Extract all hidden input fields (id, export, confirm, uuid)
        # This regex finds <input type="hidden" name="KEY" value="VALUE">
        data = {}
        hidden_inputs = re.findall(r'name="([^"]+)"\s+value="([^"]+)"', res.text)
        for key, value in hidden_inputs:
            data[key] = value

        if "confirm" in data:
            # Step 3: Use the NEW endpoint found in your HTML snippet
            confirm_url = "https://drive.usercontent.google.com/download"
            res = session.get(confirm_url, params=data, stream=True)
        else:
            raise RuntimeError("Check file permissions: 'Anyone with link' must be enabled.")

    # Step 4: Stream the actual file to disk
    res.raise_for_status()
    print(f"Downloading to {dest_path}...")
    
    with open(dest_path, "wb") as f:
        for chunk in res.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    
    print("Download complete.")
    return dest_path

# Usage
#download_google_drive_file(
#    "https://drive.google.com/file/d/1NLlVKTpWDo92NgnFapm1rk-N9yB_-6u3/view?usp=share_link",
#    "OUTCAR_collect"
#)
import sys

if len(sys.argv) < 2:
    raise SystemExit("Usage: python download.py <google-drive-url>")

download_google_drive_file(sys.argv[1], "OUTCAR_collect")

