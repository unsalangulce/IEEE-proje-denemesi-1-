from pathlib import Path

def check_h5ad_format(file_path, max_size_mb=100):
    
    if not file_path.is_file():
        print(f"File not found: {file_path}")
        return False
    
    if file_path.suffix != '.h5ad':
        print(f"The file {file_path} is not in .h5ad format.")
        return False
    
    file_size_mb = file_path.stat().st_size / (1024 * 1024)  # Convert bytes to megabytes
    if file_size_mb > max_size_mb:
        print(f"The file {file_path} is too big to be accepted. Size: {file_size_mb:.2f} MB")
        return False
    
    print(f"The file {file_path} is in .h5ad format and its size is acceptable. Size: {file_size_mb:.2f} MB")
    return True
