import requests
import zipfile
import io
import os

def main():
    """
    Downloads the tutorial notebooks from the SCOPE GitHub repository (dev branch)
    and extracts them into a local 'Scope_tutorials' folder.
    """

    repo_url = "https://github.com/QTC-IQAC/Scope/archive/refs/heads/dev.zip"
    target_dir = os.path.abspath("Scope_tutorials")

    print(f"Downloading tutorials from:\n  {repo_url}")

    r = requests.get(repo_url)
    if r.status_code != 200:
        raise RuntimeError(f"Failed to download archive (status {r.status_code})")

    # Unzip in memory
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        members = [m for m in z.namelist() if "tutorials/" in m]
        if not members:
            raise RuntimeError("No tutorials found in the downloaded archive!")
        os.makedirs(target_dir, exist_ok=True)
        for m in members:
            if m.endswith("/"):
                continue
            filename = os.path.basename(m)
            if filename:  # avoid directory entries
                with z.open(m) as source, open(os.path.join(target_dir, filename), "wb") as target:
                    target.write(source.read())

    print(f"Tutorials downloaded successfully to:\n  {target_dir}")