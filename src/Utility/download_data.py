import os
import requests

def download_file(url, output_path):
    """Download a file from the internet."""
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(1024):
                f.write(chunk)
        print(f"Downloaded: {output_path}")
    else:
        print(f"Failed to download: {url}")

# Example usage:
if __name__ == "__main__":
    url = "https://example.com/reference_genome.fasta.gz"
    output_path = "/path/to/downloaded/file.fasta.gz"
    download_file(url, output_path)
